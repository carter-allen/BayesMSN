#' Compute WAIC for fitted MSN model
#'
#' This computes WAIC for fitted MSN models to be used in model selection. 
#' @param fit object returned by fit_msn()
#' @keywords Bayesian, MSN, mixture model, clustering
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples
#' data(example1_data)
#' fit1 <- fit_msn(Y = example1_data$Y,
#'                 X = example1_data$X,
#'                 W = example1_data$W,
#'                 K = 3,
#'                 nsim = 100,
#'                 burn = 0)
#' waic(fit1)


waic <- function(fit)
{
    BETA.list = fit$BETA
    DELTA = fit$DELTA
    SIGMA.list = fit$SIGMA
    PSI.list = fit$PSI
    Z = fit$Z
    Y = fit$Y
    X = fit$X
    W = fit$W
    T = fit$T
    
    n = nrow(Y)
    K = length(BETA.list)
    J = ncol(Y)
    p = ncol(X)
    r = ncol(W)
    S = nrow(Z)
    
    lppd <- 0
    pwaic <- 0
    
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    for(i in 1:n)
    {
        y_i <- Y[i,]
        x_i <- X[i,]
        z_i <- Z[,i]
        p_yi <- rep(0,S)
        for(s in 1:S)
        {
            k <- z_i[s]
            t_i <- T[s,i]
            
            beta_k <- BETA.list[[k]]
            beta_k_s <- beta_k[s,]
            beta_ks <- matrix(beta_k_s,
                              nrow = p,
                              ncol = J,
                              byrow = FALSE)
            
            sigma_k <- SIGMA.list[[k]]
            sigma_k_s <- sigma_k[s,]
            sigma_ks <- matrix(sigma_k_s,
                               nrow = J,
                               ncol = J,
                               byrow = TRUE)
            
            psi_k <- PSI.list[[k]]
            psi_ks <- psi_k[s,]
            
            zeta_iks <- x_i %*% beta_ks
            
            mu_iks <- zeta_iks + t_i %*% psi_ks
            
            p_yi[s] <- mvtnorm::dmvnorm(y_i,mu_iks,sigma_ks)
        }
        setTxtProgressBar(pb, i)
        
        lppd <- lppd + log(mean(p_yi, na.rm = TRUE))
        pwaic <- pwaic + var(log(p_yi), na.rm = TRUE)
    }
    close(pb)
    WAIC <- -2*(lppd - pwaic)
    
    return(WAIC)
}