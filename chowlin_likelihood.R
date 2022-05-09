chowlin_likelihood <- function(Y,X,vcov) {
  
  n_l = dim(Y)[1]
  
  # Simplification and Cholesky factorization of the Sigma 
  
  Uchol <- chol(vcov)
  Lchol <- t(Uchol)
  
  # Preconditioning the variables
  
  X_F <- solve(Lchol) %*% X
  Y_F <- solve(Lchol) %*% Y
  
  # Estimate betaHat_0 using GLS assuming Sigma with rho
  
  betaHat <- solve(t(X_F) %*% X_F) %*% t(X_F) %*% Y_F 
  
  # Obtain the residuals using betaHat_0
  
  u_l_sim <- Y - X %*% betaHat
  
  # Preconditioning for the LF function
  
  u_l_sim_F <- solve(Lchol) %*% u_l_sim
  
  # Calculate the likelihood function
  
  LF <- -(n_l/2)*log(2*pi)-1/2*log(det(Lchol %*% Uchol)) - n_l/2 - n_l/2*log((t(u_l_sim_F) %*% u_l_sim_F)/n_l)
  
  return(LF)
  
}