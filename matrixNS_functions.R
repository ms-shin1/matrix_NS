Covmat <- list(
  "hub"=function(p,rho){
    Sig1 <- diag(p)
    kmax=round(p/10)
    for(k in 1:kmax){
      i <- 10*(k-1)+1
      jmin <- 10*(k-1)+2
      jmax <- 10*k
      for(j in jmin:jmax){
        Sig1[i,j] <- rho
        Sig1[j,i] <- rho
      }
    }
    eigenmin = abs(min(eigen(Sig1)$values))
    Sig = Sig1+ (eigenmin+0.05)*diag(p)
    return(Sig)
  },
  "band"=function(p,rho){
    Sig <- diag(p)
    for(i in 1:(p-1)){
      Sig[i,i+1] <- rho
      Sig[i+1,i] <- rho
    }
    for(i in 1:(p-2)){
      Sig[i,i+2] <- rho/2
      Sig[i+2,i] <- rho/2
    }
    return(Sig)
  },
  "random"=function(p,rho1,rho2){
    Sig1 <- diag(p)
    prob <- min(0.05,5/p)
    u <- runif(p*(p-1)/2,min=rho1, max=rho2)
    delta <- rbinom(p*(p-1)/2, size=1, prob)
    index <- 1
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        Sig1[i,j] <- (u*delta)[index]
        Sig1[j,i] <- (u*delta)[index]
        index <- index + 1
      }
    }
    eigenmin = abs(min(eigen(Sig1)$values))
    Sig = Sig1+ (eigenmin+0.05)*diag(p)
    return(Sig)
  }
) 

SparsityLev = function(A){
  nonzerocount <- function(x) {
  return (sum(abs(x) > 10^(-10)))
  }
  s1 <- apply(A, 1, nonzerocount)
  return(max(s1))
}

matrixrowsel <- function(X,a,lam){
  lassofit <- glmnet(y = X[,a],x = X[,-a], alpha = 1, lambda = lam)
  lassocoef <- coef(lassofit)
  nonzero <- which(lassocoef[-1] != 0)
  return(nonzero)
}

matrixrowNS <- function(X,lam){
  a<-1:nrow(X)
  sapply(a,FUN=matrixrowsel,X=X,lam=lam)
}

matrixrowsel(mtcars,1,1)
X<-mtcars
a<-1
lam<-1
a<-1:ncol(X)

