Beta <- function(y,x,ndraw){
  n <- nrow(x)
  K <- ncol(x)
  Beta <- matrix(c(0), K, ndraw)
  for(m in 1:ndraw){
    V <- rexp(n,1)
    sV <- sqrt(V)
    Ytilde <- sV * y
    Xtilde <- sV * x
    Beta[,m] <- solve((t(Xtilde)%*% Xtilde)) %*% (t(Xtilde) %*% Ytilde)
  }
  return(Beta)
}

