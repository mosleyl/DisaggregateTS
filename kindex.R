
k.index <- function(matrix, n_l) {
  
  count <- apply(matrix, 1, function(x) {sum(x != 0)})
  if(max(count) > n_l/2) {
    kindex <- min(which(count > n_l/2))
  }
  else {
    kindex <- nrow(matrix)
  }
  
  return(kindex)
}