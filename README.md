BayesMSN
===

This package provides implementations of the Bayesian Multivariate Skew-Normal (MSN) mixture models presented in Allen et al. (2020). The Bayesian MSN model is designed to uncover latent clusters in multivariate/longitudinal data that may exibit skewness, e.g., repeated measures infant development scores. 

Installation
===
```
library(devtools)
install_github('carter-allen/BayesMSN')
```

Usage
===
```
library(BayesMSN)
data(example1_data)
fit1 <- fit_msn(Y = example1_data$Y,
                X = example1_data$X,
                W = example1_data$W,
                K = 3,
                nsim = 10,
                burn = 0)
```