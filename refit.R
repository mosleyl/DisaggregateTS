refit <- function(X, Y, beta) {
  
  p = ncol(X)
  active <- which(beta != 0 )
  X_active <- X[,active]
  lm_fit <- lm(Y ~ 0 + X_active)
  betahat_ols <- lm_fit$coefficients
  betahat_debias <- rep(0,p)
  betahat_debias[active] <- betahat_ols
  
  return(betahat_debias)
  
}