#' Fit MSN mixture model 
#'
#' This function fits the MSN mixture model without imputation
#' @param Y matrix of non-missing outcome data
#' @param X matrix of predictors used in the MSN model component
#' @param W matrix of predictors used in the multinomial logit clustering component
#' @param K number of clusters to fit. Use WAIC to choose best K.
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
#' fit1 <- fit_msn(Y = example1_data$Y,
#'                 X = example1_data$X,
#'                 W = example1_data$W,
#'                 K = 3,
#'                 nsim = 10,
#'                 burn = 0)

fit_msn <- function(Y,X,W,K,nsim,burn)
{
    t = truncnorm::rtruncnorm(nrow(X),a = 0) # initialize the random truncnorm latent vars
    Xstar <- cbind(X,t) # combine with design matrix
    
    # Define parameters
    h = K # number of clusters
    k = ncol(Y) # number of repeated measures outcomes
    n = nrow(Y)  # number of subjects
    p = ncol(X) # number of MSN regression covariates
    v = ncol(W) # number of multinomial reg. covariates
    
    # Define priors 
    V0 = diag(1,k) # used in matnorm updates
    V0.list = vector("list",h)
    V0.list[1:h] = list(V0)
    L0.list = vector("list",h) # used in matnorm updates
    L0.list[1:h] = list(diag(1,p+1))
    B0.list = vector("list",h) # prior mean for MSN reg. coeffs
    B0.list[1:h] = list(matrix(0,nrow = p+1,ncol = k))
    delta0 = rep(0,v) # prior mean for delta coefficients (multinomial regression)
    S0 = diag(1,v) # prior covariance for delta coefficients (multinomial regression)
    nu0 = rep(6,h) 
    
    # Define inits
    pi = rep(1/h,h) # used to initialize z (uniformly)
    z = sample(1:h,n,replace = TRUE,prob = pi) # initialize z
    Bn.list = B0.list # make storage for B
    Ln.list = L0.list # make storage for L
    nun = nu0 # used in updates of Sigma
    Vn.list = V0.list # used in updates of Sigma
    beta.list = vector("list",h) #  storage for betas
    beta.mle = solve(t(X) %*% X) %*% t(X) %*% Y # mle of betas
    beta.list[1:h] = list(beta.mle) # inits as mles
    sig2.list = vector("list",h) # storage for sigmas
    sig2.mle = mlest(Y - X %*% beta.mle)$sigmahat # mles of sigmas
    sig2.list[1:h] = list(sig2.mle) # inits as mles
    omega.list = V0.list # storage for omegas
    psi.list = vector("list",h) # init for psis
    psi.list[1:h] = list(rep(0,k))
    alpha.list = psi.list  # init for alphas 
    delta = matrix(0,nrow = v, ncol = h-1) # storage for betas

    n.iter <- nsim - burn # number of saved iterations
    Z = matrix(0,nrow = n.iter,ncol = n) # large matrix where each row is the value of z at a specific iteration
    T = matrix(0,nrow = n.iter,ncol = n) # storage for t_i's
    PI = matrix(0,nrow = n.iter,ncol = h) # matrix w/ each row as pi vector at each iteration
    BETA.list = vector("list",h) # storage for Beta (vectorized by row)
    SIGMA.list = vector("list",h) # storage for S matrix (vectorized by row)
    PSI.list = vector("list",h) # storage for psi vector
    DELTA = matrix(0,nrow = n.iter,ncol = v*(h-1)) # storage for delta matrix
    OMEGA.list = vector("list",h) # storage for omega matrix (vectorized by row)
    ALPHA.list = vector("list",h) # storage for alphas
    
    # MCMC
    start.time<-proc.time()
    print(paste("Started MCMC of",nsim,"iterations with a burn in of",burn))
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    for(i in 1:nsim)
    {
        # Step 1:
        # Update pi1,...,piK 
        eta <- cbind(rep(0,n),W%*%delta)
        pi <- exp(eta)/(1+apply(as.matrix(exp(eta[,-1])),1,sum))
        
        # Step 2:
        # Update z_j from multinomial conditional posterior
        # for each observed y, compute the probability of belonging in each class
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
                p.j[l] <- dmvnorm(y.j,mu.j,sig2.l) # get mvn density
            }
            pi.j <- pi[j,]
            p.z <- pi.j*p.j/sum(pi.j*p.j)
            z[j] <- which(rmultinom(1,1,p.z) == 1)
        }
        # convert z to factor so that empty classes are included as 0 counts below and not dropped
        z <- factor(z,levels = 1:h)
        # # compute the sample size within each class
        # # n should be a vector of length h
        # # if there is an empty class then the program should stop
        n.z <- as.vector(unname(table(z))) # gives the number of members currently in each class
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
        
        # Step 3:
        # Update multinomial regression parameters
        W <- as.matrix(W) # enforce W is a matrix
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
        
        # Step 4: 
        # Update class-specific regression parameters
        T.i <- rep(0,n)
        for(l in 1:h) # loop through each cluster
        {
            X.l <- as.matrix(X[z == l,]) # all covariates for those in class l
            Y.l <- Y[z == l,] # all outcomes for those in class l
            n.l <- sum(z == l) # current class specific sample size
            psi.l <- psi.list[[l]] # current class specific psis 
            B0.l <- B0.list[[l]] # priors B0
            V0.l <- V0.list[[l]] # priors V0
            L0.l <- L0.list[[l]] # priors L0
            nu0.l <- nu0[l] # priors nu0
            nun.l <- nun[l] # current nu
            Bn.l <- Bn.list[[l]] # current B
            Vn.l <- Vn.list[[l]] # current V
            Ln.l <- Ln.list[[l]] # current L 
            sig2.l <- sig2.list[[l]] # current sigma matric
            beta.l <- beta.list[[l]] # current beta matrix
            betastar.l <- rbind(beta.l,psi.l) # current betastar
            
            # Update t
            A <- solve(1 + t(psi.l) %*% solve(sig2.l) %*% psi.l)
            t.l <- rep(0,n.l)
            for(j in 1:n.l)
            {
                ai <- A %*% t(psi.l) %*% solve(sig2.l) %*% t((t(Y.l[j,]) - t(X.l[j,]) %*% beta.l))
                t.l[j] <- truncnorm::rtruncnorm(n = 1,a = 0, b = Inf,mean = ai, sd = sqrt(A))
            }
            Xstar.l <- cbind(X.l,t.l) # add updated t's back to Xstar matrix
            T.i[z == l] <- t.l # save T's 
            
            # Update sigma
            # Same as matrix normal regression update
            Vn.l <- V0.l + t(Y.l - Xstar.l %*% Bn.l) %*% (Y.l - Xstar.l %*% Bn.l) + t(Bn.l - B0.l) %*% L0.l %*% (Bn.l - B0.l)
            nun.l <- nu0.l + n.l
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