#' Fit MST mixture model 
#'
#' This function fits the MST mixture model without imputation
#' @param Y matrix of non-missing outcome data
#' @param X matrix of predictors used in the MSN model component
#' @param W matrix of predictors used in the multinomial logit clustering component
#' @param K number of clusters to fit. Use WAIC to choose best K.
#' @param nu degrees of freedom parameter for skew-t model.
#' @param nsim number of MCMC iterations to perform
#' @param burn number of initial MCMC iterations to discard. The function will return nsim - burn total posterior samples.
#' @keywords Bayesian, MSN, mixture model, clustering
#' @export
#' @import sn
#' @import truncnorm
#' @import mvtnorm 
#' @import mvnmle
#' @import Hmisc 
#' @import coda
#' @importFrom MCMCpack riwish
#' @import matrixsampling
#' @import nnet
#' @importFrom BayesLogit rpg
#' @importFrom stats rmultinom
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(example1_data)
#' fit1 <- fit_mst(Y = example1_data$Y,
#'                 X = example1_data$X,
#'                 W = example1_data$W,
#'                 K = 3,
#'                 nsim = 10,
#'                 burn = 0)

fit_mst <- function(Y,X,W,K,
                    nu = 3,
                    nsim = 1000,
                    burn = 0)
{
  t.init <- rtruncnorm(nrow(X),0,Inf)
  Xstar <- cbind(X,t.init)
  Y = as.matrix(Y)
  Xstar = as.matrix(Xstar)
  W = as.matrix(W)
  X = Xstar[,1:(ncol(Xstar)-1)]
  
  # Define parameters
  h = K
  k = ncol(Y)
  n = nrow(Y)
  p = ncol(X)
  v = ncol(W)
  
  # Define priors 
  V0 = diag(1,k)
  V0.list = vector("list",h)
  V0.list[1:h] = list(V0)
  L0.list = vector("list",h)
  L0.list[1:h] = list(diag(1,p+1))
  B0.list = vector("list",h)
  B0.list[1:h] = list(matrix(0,nrow = p+1,ncol = k))
  delta0 = rep(0,v) # prior mean for delta coefficients (multinomial regression)
  S0 = diag(1,v) # prior covariance for delta coefficients (multinomial regression)
  nu0 = rep(6,h) 
  
  Y.true <- Y
  if(any(is.na(Y)))
  {
    
    Y <- mice::complete(mice::mice(as.data.frame(Y), m = 1))
    Y <- as.matrix(Y)
  }
  na.mat <- is.na(Y.true)
  
  # Define inits
  pi = rep(1/h,h)
  print("Randomly initializing cluster allocations")
  z = sample(1:h,n,replace = TRUE,prob = pi)
  Bn.list = B0.list
  Ln.list = L0.list
  nun = nu0
  Vn.list = V0.list
  beta.list = vector("list",h)
  beta.mle = solve(t(X) %*% X) %*% t(X) %*% Y
  beta.list[1:h] = list(beta.mle)
  sig2.list = vector("list",h)
  sig2.mle = mlest(Y - X %*% beta.mle)$sigmahat
  sig2.list[1:h] = list(sig2.mle)
  omega.list = V0.list
  psi.list = vector("list",h)
  psi.list[1:h] = list(rep(0,k))
  alpha.list = psi.list
  delta = matrix(0,nrow = v, ncol = h-1)
  
  n.iter <- nsim - burn # number of saved iterations
  Z = matrix(0,nrow = n.iter,ncol = n) # large matrix where each row is the value of z at a specific iteration
  T = matrix(0,nrow = n.iter,ncol = n) # storage for t_i's
  PI = matrix(0,nrow = n.iter,ncol = h) # matrix w/ each row as pi vector at each iteration
  BETA.list = vector("list",h) # storage for Beta (vectorized by row)
  SIGMA.list = vector("list",h) # storage for S matrix (vectorized by row)
  PSI.list = vector("list",h) # storage for psi vector
  DELTA = matrix(0,nrow = n.iter,ncol = v*(h-1))
  OMEGA.list = vector("list",h) # storage for omega matrix (vectorized by row)
  ALPHA.list = vector("list",h)
  Z_bank = z
  
  # MCMC
  start.time<-proc.time()
  print(paste("Started MCMC of",nsim,"iterations with a burn in of",burn))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(i in 1:nsim)
  {
    # Step 0:
    # perform conditional mv-normal imputation
    for(j in 1:n)
    {
      # perform conditional normal imputation
      na.ind <- na.mat[j,]
      y.j <- Y[j,] # outcome for observation j
      if(any(na.ind))
      {
        if(all(na.ind))
        {
          y.imp <- colMeans(Y,na.rm = TRUE)
        }
        else
        {
          a.j <- y.j[!na.ind] # observed components of y.j
          xstar.j <- Xstar[j,] # covariates for observation j
          z.j <- z[j] # current cluster id of y.j
          A.j <- sig2.list[[z.j]] # covmatrix of y.j
          betastar.j <- rbind(beta.list[[z.j]],psi.list[[z.j]]) # betstar for y.j
          mu.j <- xstar.j %*% betastar.j # mean of observation j if it belonged to each class
          
          mu.1 <- mu.j[na.ind]
          mu.2 <- mu.j[!na.ind]
          
          sig.11 <- A.j[na.ind,na.ind]
          sig.12 <- A.j[na.ind,!na.ind]
          sig.21 <- A.j[!na.ind,na.ind]
          sig.22 <- A.j[!na.ind,!na.ind]
          
          mu.cond <- mu.1 + sig.12 %*% solve(sig.22) %*% (a.j - mu.2)
          sig.cond <- sig.11 - sig.12 %*% solve(sig.22) %*% sig.21
          if(!isSymmetric(sig.cond))
          {
            sig.cond <- round(sig.cond,digits = 8)
          }
          
          y.imp <- rmvnorm(1,mu.cond,sig.cond)
        }
        y.j[na.ind] <- y.imp
        Y[j,] <- y.j
      }
    }
    Y <- as.matrix(Y)
    
    # Step 3:
    # Update pi1,...,piK 
    eta <- cbind(rep(0,n),W%*%delta)
    pi <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
    # End Step 3
    
    # Step 1A:
    # Update z_j from multinomial conditional posterior
    # for each observed y, compute the probability of belonging in each class
    # this is done iteratively in a loop but there is probably a better way
    for(j in 1:n)
    {
      y.j <- Y[j,] # outcome for observation j
      xstar.j <- Xstar[j,] # covariates for observation j
      p.j <- rep(0,h)
      for(l in 1:h)
      {
        betastar.l <- rbind(beta.list[[l]],psi.list[[l]])
        mu.j <- xstar.j %*% betastar.l # mean of observation j if it belonged to each class
        sig2.l <- sig2.list[[l]]
        p.j[l] <- dmvnorm(y.j,mu.j,sig2.l)
      }
      pi.j <- pi[j,]
      p.z <- pi.j*p.j/sum(pi.j*p.j)
      # print(p.z)
      z[j] <- which(rmultinom(1,1,p.z) == 1) # this is very slow 
      # z <- c # for testing purposes hold z to true value
      # convert z to factor so that empty classes are included as 0 counts below and not dropped
      z = factor(z,levels = 1:h)
      # Z_bank = rbind(Z_bank,z)
    }
    # # compute the sample size within each class
    # # n should be a vector of length h
    # # if there is an empty class then the program should stop
    n.z <- as.vector(unname(table(z))) # gives the number of members currently in each class
    # print(n.z/n)
    if(any(n.z == 0)) # check for empty clusters, if so stop
    {
      print("MCMC terminated due to an empty class.")
      break
    }
    if(any(n.z == 1)) # check for singleton clusters, if so stop
    {
      print("MCMC terminated due to a singleton class.")
      break
    }
    # End Step 1A
    
    # Step 1B:
    # Update multinomial regression parameters
    W <- as.matrix(W)
    for(l in 1:(h-1))
    {
      delta.l <- delta[,l]
      delta.notl <- delta[,-l]
      u.l <- 1*(z == l+1)
      c.l <- log(1 + exp(rowSums(W %*% delta.notl)))
      eta <- W %*% delta.l - c.l
      w <- rpg(n, 1, eta)
      z.l <- (u.l - 1/2)/w + c.l 
      
      V <- solve(S0 + crossprod(W*sqrt(w)))  
      M <- V %*% (S0 %*% delta0 + t(w*W) %*% z.l)
      delta.l <- c(rmvnorm(1,M,V))
      delta[,l] <- delta.l
    }
    # End Step 1B
    
    # Step 2: 
    # Update class-specific regression parameters
    T.i <- rep(0,n)
    for(l in 1:h) # loop through each cluster
    {
      X.l <- as.matrix(X[z == l,]) # all covariates for those in class l
      Y.l <- Y[z == l,] # all outcomes for those in class l
      n.l <- sum(z == l) # current class specific sample size
      psi.l <- psi.list[[l]] # current class specific psis 
      B0.l <- B0.list[[l]]
      V0.l <- V0.list[[l]]
      L0.l <- L0.list[[l]]
      nu0.l <- nu0[l]
      nun.l <- nun[l]
      Bn.l <- Bn.list[[l]]
      Vn.l <- Vn.list[[l]]
      Ln.l <- Ln.list[[l]]
      sig2.l <- sig2.list[[l]]
      beta.l <- beta.list[[l]]
      betastar.l <- rbind(beta.l,psi.l)
      
      # Update t
      A <- solve(1 + t(psi.l) %*% solve(sig2.l) %*% psi.l)
      t.l <- rep(0,n.l)
      d.l <- rep(0,n.l) # empty latend gamma
      Xstar.l <- cbind(X.l,t.l) # add updated t's back to Xstar matrix
      for(j in 1:n.l)
      {
        ai <- A %*% t(psi.l) %*% solve(sig2.l) %*% t((t(Y.l[j,]) - t(X.l[j,]) %*% beta.l))
        t.l[j] <- rtruncnorm(n = 1,a = 0, b = Inf,mean = ai, sd = sqrt(A))
        Xstar.l[j,p+1] <- t.l[j]
        ej <- Y.l[j,] - Xstar.l[j,] %*% betastar.l
        d.l[j] <- rgamma(n = 1,(nu+k+1)/2,(nu+t.l[j]^2+ej %*% solve(sig2.l) %*% t(ej))/2)
      }
      Y.l <- d.l * Y.l
      Xstar.l <- d.l * Xstar.l
      T.i[z == l] <- t.l
      
      # Update sigma
      # Same as matrix normal regression update
      Vn.l <- V0.l + t((Y.l) - Xstar.l %*% Bn.l) %*% (Y.l - Xstar.l %*% Bn.l) + t(Bn.l - B0.l) %*% L0.l %*% (Bn.l - B0.l)
      nun.l <- nu0.l + n.l + p + 1
      sig2.l <- riwish(nun.l,Vn.l)
      sig2.list[[l]] <- sig2.l
      omega.l <- sig2.l + outer(psi.l,psi.l)
      omega.list[[l]] <- omega.l
      
      # Update beta
      # Same as matrix normal regression update with added psi
      Bn.l <- solve(t(Xstar.l) %*% Xstar.l + L0.l) %*% (t(Xstar.l) %*% Y.l + L0.l %*% B0.l)
      Bn.list[[l]] <- Bn.l
      Ln.l <- t(Xstar.l) %*% Xstar.l + L0.l
      Ln.list[[l]] <- Ln.l
      betastar.l <- rmatrixnormal(1,Bn.l,solve(Ln.l), sig2.l, checkSymmetry = FALSE)
      betastar.l <- matrix(betastar.l,nrow = p+1, ncol = k)
      beta.l <- betastar.l[1:p,]
      psi.l <- betastar.l[p+1,]
      beta.list[[l]] <- beta.l
      psi.list[[l]] <- psi.l
    }
    # End step 2
    
    
    # Store results
    if (i > burn)
    {
      j <- i-burn
      Z[j,] <- z
      PI[j,] <- table(z)/n
      DELTA[j,] <- c(delta)
      T[j,] <- T.i
      for(l in 1:h)
      {
        # ALPHA.list[[l]] <- rbind(ALPHA.list[[l]],alpha.list[[l]])
        OMEGA.list[[l]] <- rbind(OMEGA.list[[l]],c(omega.list[[l]]))
        BETA.list[[l]] <- rbind(BETA.list[[l]],c(beta.list[[l]]))
        SIGMA.list[[l]] <- rbind(SIGMA.list[[l]],c(sig2.list[[l]]))
        PSI.list[[l]] <- rbind(PSI.list[[l]],psi.list[[l]])
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  run.time<-proc.time()-start.time
  print(paste("Finished MCMC after",run.time[1],"seconds"))
  ret_list = list(BETA = BETA.list,
                  DELTA = DELTA,
                  SIGMA = SIGMA.list,
                  PSI = PSI.list,
                  T = T,
                  Z = Z,
                  Y = Y,
                  X = X,
                  W = W)
  return(ret_list)
}