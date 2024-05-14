rm(list=ls())
library(jointMeanCov)
library(matrixNormal)

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

set.seed(1234)
X<-rmatnorm(s=1,zero,U,V)
out <- GeminiB(X, rowpen=0.01)
out$B.hat


U <- Covmat[['hub']](p,rho)
V <- Covmat[['band']](q,rho)
R <- Covmat[['random']](200,0.3,0.7)

SparsityLev(solve(V))

solve(U)
solve(V)
solve(R)

memory.size(max = TRUE) # OS로부터 R이 사용 가능한 메모리 
memory.size(max = FALSE) # 현재 사용중인 메모리 크기 
memory.limit(size = NA) # 컴퓨터의 최대 메모리 한계치 
memory.limit(size = 50000) # 컴퓨터의 최대 메모리 한계치 약 49GB로 높이기
