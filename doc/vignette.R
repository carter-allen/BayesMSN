## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)

## ---- eval = FALSE------------------------------------------------------------
#  devtools::install_github("carter-allen/BayesMSN")

## -----------------------------------------------------------------------------
library(BayesMSN)

## -----------------------------------------------------------------------------
str(example1_data)

## ---- results='hide'----------------------------------------------------------
set.seed(1801)
fit3 <- fit_msn(Y = example1_data$Y,
                X = example1_data$X,
                W = example1_data$W,
                K = 3,
                nsim = 100,
                burn = 0)

## ---- results='hide'----------------------------------------------------------
fit3_waic <- waic(fit3)

## -----------------------------------------------------------------------------
fit3_waic

## -----------------------------------------------------------------------------
col_summarize(fit3$BETA[[1]])

## -----------------------------------------------------------------------------
col_summarize(fit3$DELTA)

## ---- results='hide'----------------------------------------------------------
fit3_mnar <- fit_msn_impute(example2_data$Y,
                     example2_data$X,
                     example2_data$W,
                     Q = example2_data$Q,
                     K = 3,
                     nsim = 100)

## -----------------------------------------------------------------------------
G1 <-  fit3_mnar$GAMMA[[1]]
G1 <- as.data.frame(G1)
colnames(G1) <- c("g11", "g12", "g13")

## -----------------------------------------------------------------------------
G1 %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>%
  ggplot(.,aes(x = variable, y = value)) + 
  geom_boxplot() + 
  ggtitle("Cluster 1 missing data model") + 
  theme_classic()

## ---- results='hide'----------------------------------------------------------
fit3_MST_3df <- fit_mst(Y = example1_data$Y,
                    X = example1_data$X,
                    W = example1_data$W,
                    K = 3,
                    nu = 3,
                    nsim = 100,
                    burn = 0)

## ----  results='hide'---------------------------------------------------------
fit3_MST_3df_waic <- waic(fit3_MST_3df)

## -----------------------------------------------------------------------------
fit3_MST_3df_waic

