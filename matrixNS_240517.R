rm(list=ls())
library(jointMeanCov)
library(matrixNormal)
library(glmnet)
library(expm)

wd = getwd()
source(paste0(wd,"/matrixNS_functions.R"))


p<-20
q<-20
n<-100
U<-diag(p)
V<-diag(q)

lambda<-0.1
rho <- 0.5

Uinv <- Covmat[['hub']](p,rho)
Vinv <- Covmat[['band']](q,rho)
U <- solve(Uinv)
V <- solve(Vinv)

set.seed(1234)
Xvec <- rnorm(p*q,mean=0,sd=1)
X <- sqrtm(U) %*% matrix(Xvec, nrow=p) %*% sqrtm(V)
matrixcolNS(X,0.6)

R <- Covmat[['random']](200,0.3,0.7)

