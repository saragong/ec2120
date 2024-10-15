beta <- function(y,x,A,Phi,weight){
  n <- nrow(y)
  M <- ncol(y)
  q <- ncol(x)
  K <- ncol(A)
  sqw <- sqrt(weight)
  xw <- sqw * x
  yw <- sqw * y
  Pi <- solve( ( t(xw) %*% xw ) ) %*% ( t(xw) %*% y )
  uw <- yw - xw %*% Pi
  Sigma <- (t(uw) %*% uw) / sum(weight)
  Sigma_inv <- solve(Sigma)
  d1 <- matrix(c(0),K,K)
  d2 <- matrix(c(0),K,1)
  for(i in 1:n){
    rwi <- (diag(M) %x% t(xw[i,])) %*% A
    d1 <- d1 + t(rwi) %*% Phi %*% rwi
    d2 <- d2 + t(rwi) %*% Phi %*% yw[i,]
  }
  beta <- solve((t(d1) %*% d1)) %*% (t(d1) %*% d2)
  return(beta)
}