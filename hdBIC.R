hdBIC <- function(X, Y, covariance, beta) {
  
  n_l <- length(Y)
  support <- sum(beta != 0)
  u <- Y-X%*%beta
  log.lik <- -(n_l-support)/2 - n_l/2*log(2*pi/(n_l-support)*crossprod(u)) - log(det(covariance))/2
  
  BIC <- -2*(log.lik) + log(n_l)*support
  #AICc <- -2*(log.lik) + 2*support + 2*support*(support+1)/(n_l-support-1)
  
  return(BIC)
  
}
