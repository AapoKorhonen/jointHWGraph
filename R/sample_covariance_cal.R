
sample_covariance_cal <- function(data_matrix){
  
  n <- nrow(data_matrix)
  
  X_centered <- scale(data_matrix, center = T, scale = F)
  
  S <- t(X_centered) %*% X_centered / (n-1)
  return(S)
}