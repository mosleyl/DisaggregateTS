ARcov_lit <- function(rho, n) {
  
  diags <- list(rep(1, times = n), rep(-rho, times = n-1))
  H_r_t <-  bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
  H_r <- t(H_r_t) 
  
  Sigma <- solve(Delta_t %*% H_r_t %*% H_r %*% Delta)
  
}