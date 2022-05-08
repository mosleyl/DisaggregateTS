disaggregate <- function(Y, X = matrix(data = rep(1, times = nrow(Y)), nrow = nrow(Y)), aggMat = 'sum', aggRatio = 4, method = 'Denton-Cholette', Denton = 'first'){
  
  
  if(is.matrix(X) == FALSE || is.matrix(Y) == FALSE){
    
    stop("X and Y must be a matrices! \n")
    
  }
  
  if(!(method == 'Denton' || method == 'Denton-Cholette' || method == 'Chow-Lin' || method == 'Fernandez' || method == 'Litterman' || method == 'spTD' || method == 'adaptive-spTD')) {
    
    stop("Wrong method inputted \n")
    
  }
  
  nl = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  nfull = aggRatio*nl
  extr = n - nfull # number of extrapolations
  
  if(nfull > n) {
    
    stop("X does not have enough observations. \n")
    
  }
  
  
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
  
  
  if(method == 'Denton-Cholette' || method == 'Denton'){
    
    # Difference between the low-frequency and the transformed high-frequency series.
    
    u_l <- Y - Xl
    
    # First difference matrix
    
    diags <- list(rep(1, times = n), rep(-1, times = n-1))
    Delta_t <-  bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
    Delta <- t(Delta_t)
    
    if(method == 'Denton'){
      
      # Absolute difference matrix
      
      if(Denton == 'abs'){
        
        Sigma <- diag(n)
        
        # First difference Sigma
        
      }else if(Denton == 'first'){
        
        Sigma <- solve(Delta_t %*% Delta)
        
        # Second difference Sigma
        
      }else if(Denton == 'second'){
        
        Sigma <- solve(Delta_t %*% Delta_t %*% Delta %*% Delta) 
        
        # Proportional difference Sigma
        
      }else if(Denton == 'prop'){
        
        Sigma <- solve(solve(Diagonal(n, X)) %*% Delta_t %*% Delta %*% solve(Diagonal(n, X)))
        
      }
      
      # The distribution matrix
      
      D <- Sigma %*% t(C) %*% solve(C %*% Sigma %*% t(C))
      
      # Generate the high-frequency series
      
      y <- X + (D %*% u_l) 
      
      # Unnecessary parameter outputs
      
      rho_opt <- NaN
      betaHat_opt <- NaN
      
    }else if(method == 'Denton-Cholette'){
      
      # Removed the first roq of the Delta matrix
      
      Delta_DC <- Delta[2: nrow(Delta), ]
      
      # The distribution matrix
      
      D <- Delta_DC
      
      # Generate the high-frequency series
      
      y <- X + (D %*% u_l) 
      
      # Unnecessary parameter outputs
      
      rho_opt <- NaN
      betaHat_opt <- NaN
      
    }
    
  }else if(method == 'Fernandez') {
    
    diags <- list(rep(1, times = n), rep(-1, times = n-1))
    Delta_t <- bandSparse(n, k = 0:1, diagonals = diags, symmetric = FALSE)
    Delta <- t(Delta_t) 
    
    # H(rho) matrix is first an nxn identity matrix
    
    Sigma_opt <- solve(Delta_t %*% Delta)
    
    # Simplification and Cholesky factorization of the Sigma_opt 
    
    Uchol_opt <- chol(C %*% Sigma_opt %*% t(C))
    Lchol_opt <- t(Uchol_opt)
    
    # Preconditioning the variables
    
    X_F_opt <- solve(Lchol_opt) %*% Xl
    Y_F_opt <- solve(Lchol_opt) %*% Y  
    
    # First estimate betaHat_opt using OLS assuming Sigma = (Delta'Delta)^{-1}
    
    betaHat_opt <- solve(t(X_F_opt) %*% X_F_opt) %*% t(X_F_opt) %*% Y_F_opt 
    
    # Obtain the residuals using betaHat_1
    
    u_l <- Y - Xl %*% betaHat_opt
    
    # The distribution matrix
    
    D <- Sigma %*% t(C) %*% solve(C %*% Sigma %*% t(C))
    
    # Generate the high-frequency series
    
    y <- X %*% betaHat_opt + (D %*% u_l)
    
    rho_opt <- NaN
    
  }else {
    
    if(method == 'Chow-Lin') {
      
      Objective <- function(rho) {
        -chowlin_likelihood(
          Y = Y, X = Xl, vcov = C %*% ARcov(rho,n) %*% t(C)
        )
      }
    }else if(method == 'Litterman') {
      
      Objective <- function(rho) {
        -chowlin_likelihood(
          Y = Y, X = Xl, vcov = C %*% ARcov_lit(rho,n) %*% t(C)
        )
      }
      
    }else if(method == 'spTD' || method == 'adaptive-spTD') {
      
      Objective <- function(rho) {
        sptd_BIC(
          Y = Y, X = Xl, vcov = C %*% ARcov(rho,n) %*% t(C)
        )
      }
      
    }
    
    # optimise for best rho 
    
    optimise_rho <- optimize(Objective,
                             lower = 0, upper = 0.999, tol = 1e-16,
                             maximum = FALSE
    )
    
    rho_opt = optimise_rho$minimum
    
    # Generate the optimal Toeplitz covariance matrix
    
    if(method == Chow-Lin){
      
      fit = chowlin(Y = Y, X = X, rho = rho_opt, aggMat = 'sum', aggRatio = 4)
      betaHat = fit$betaHat
      y = fit$y
      
    }
      
  }
  
  

}