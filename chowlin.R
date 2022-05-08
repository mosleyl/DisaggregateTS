chowlin <- function(Y, X, rho, aggMat, aggRatio, litterman = FALSE) {
  
  nl = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  nfull = aggRatio*nl
  extr = n - nfull # number of extrapolations
  
  
  # Generate the aggregation matrix C
  
  if(aggMat == 'sum'){
    
    C <- kronecker(diag(n_l), matrix(data = 1, nrow = 1, ncol = m))
    C <- cbind(C, matrix(0L, nl, extr))
    
  }else if(aggMat == 'avg'){
    
    C <- kronecker(diag(n_l), matrix(data = n_l/n, nrow = 1, ncol = m))
    C <- cbind(C, matrix(0L, nl, extr))
    
  }else if(aggMat == 'first'){
    
    C <- kronecker(diag(n_l), matrix(data = c(1, rep(0, times = m-1)), nrow = 1, ncol = m))
    C <- cbind(C, matrix(0L, nl, extr))
    
  }else if(aggMat == 'last'){
    
    C <- kronecker(diag(n_l), matrix(data = c(rep(0, times = m-1), 1), nrow = 1, ncol = m))
    C <- cbind(C, matrix(0L, nl, extr))
    
  }
  
  Xl = C %*% X
  
  if(litterman) {
    vcov = ARcov_lit(rho, n)
  }else {
    vcov = ARcov(rho, n)
  }
  
  # Simplification and Cholesky factorization of the Sigma 
  
  vcov_agg = C %*% vcov %*% t(C)
  Uchol <- chol(vcov_agg)
  Lchol <- t(Uchol)
  
  # Preconditioning the variables
  
  X_F <- solve(Lchol) %*% X
  Y_F <- solve(Lchol) %*% Y
  
  # Estimate betaHat_0 using GLS assuming Sigma with rho
  
  betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F 

  
  # The distribution matrix
  
  D <- vcov %*% t(C) %*% solve(vcov_agg)
  
  # Obtain the residuals using betaHat_1
  
  u_l <- Y - C %*% X %*% betaHat  
  
  # Generate the high-frequency series
  
  y <- X %*% betaHat + (D %*% u_l)
  
  output = list('y' = y, betaHat = 'betaHat', 'u_l' = u_l)
  

  return(output)
  
}