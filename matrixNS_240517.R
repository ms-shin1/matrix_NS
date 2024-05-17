rm(list=ls())
library(jointMeanCov)
library(matrixNormal)
library(glmnet)
library(expm)

wd = getwd()
source(paste0(wd,"/matrixNS_functions.R"))


p<-200
q<-200
n<-100
U<-diag(p)
V<-diag(q)
zero<-matrix(0,nrow=p,ncol=q)

lambda<-0.1
rho <- 0.5

Uinv <- Covmat[['hub']](p,rho)
Vinv <- Covmat[['band']](q,rho)
U <- solve(Uinv)
V <- solve(Vinv)

set.seed(1234)
Xvec <- rnorm(p*q,mean=0,sd=1)
X <- sqrtm(U) %*% matrix(Xvec, nrow=p) %*% sqrtm(V)

R <- Covmat[['random']](200,0.3,0.7)

solve(U)
solve(V)
solve(R)

memory.size(max = TRUE) # OS로부터 R이 사용 가능한 메모리 
memory.size(max = FALSE) # 현재 사용중인 메모리 크기 
memory.limit(size = NA) # 컴퓨터의 최대 메모리 한계치 
memory.limit(size = 50000) # 컴퓨터의 최대 메모리 한계치 약 49GB로 높이기
