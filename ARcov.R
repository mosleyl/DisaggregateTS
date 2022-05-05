ARcov <- function(rho,n) {
  
  sqnc <- rho^seq(0, n, by = 1)
  Omega <- toeplitz(sqnc[1: n])
  sig <- (1/(1-rho^2)) * Omega
  return(sig)
  
}