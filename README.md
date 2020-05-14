BayesMSN
===

This package provides implementations of the Bayesian Multivariate Skew-Normal (MSN) mixture models presented in Allen et al. (2020). The Bayesian MSN model is designed to uncover latent clusters in multivariate/longitudinal data that may exibit skewness, e.g., repeated measures infant development scores. 

Installation
===
```
library(devtools)
install_github('carter-allen/BayesMSN')
```

Overview
===

`BayesMSN` is an `R` package for fitting Bayesian multivariate skew-normal mixture models to longitudinal/repeated measures data that may possibly feature latent sub-clusters of longitudinal outcomes. The model is presented in full detail in Allen et al. (2020). Some key features of this model are as follows:

- Uses finite mixture modeling to explain heterogeneity among repeated measures in terms of a parsimonious set of mixture components (i.e., clusters). This feature of the model is useful if interest lies in uncovering latent clusters in longitudinal outcomes that may not be apparent from marginal trajectories.

- Accounts for the possibility of skewness or heavy tails in repeated measures by utilizing skew-normal or skew-$t$ distributions instead of traditional Gaussian distributions. 

- Allows for a wide range of possible covariance patterns among the repeated measures by utlizing unstructured covariance matrices for the multivariate skew-normal (MSN) and multivariate skew-$t$ (MST) distributions.

- Models the cluster indicators themselves using P\'olya--Gamma data augmentation to explain cluster membership as a function of covariates of interest. This is a key advantage over other clustering methods that allows for deeper understanding of cluster profiles and readily available practical interpretations of clusters. 

- Efficiently imputes intermittent missing repeated measures without globar missing at random assumptions (MAR) through the use of the cluster indicators as a discrete shared parameter to account for unobserved association between the missing mechanism and the missing values themselves. 

In summary, the `BayesMSN` package is useful for the longitudinal or multivariate data analysts who seeks a flexible model for uncovering latent sub-clusters among the responses and the ability to explain cluster membership in terms of other practically relevant data. This vignette presents three data analysis examples that showcase some of the key features of `BayesMSN`. 

___Note:___ The core functions in `BayesMSN` use Gibb's sampling to obtain samples from the posterior distributions of all model parameters. Like all Bayesian MCMC methods, care must be taken to assess convergence of the MCMC sampler. To perform these diagnostics, we suggest external packages such as `coda`, `bayesplot`, and `label.switching`. 

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

