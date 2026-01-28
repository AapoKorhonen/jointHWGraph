#' Function for thresholding the partial correlation matrix.
#'
#' @param mat The partial correlation matrix to be thresholded.
#' @param thres The threshold
#'
#' @return Returns the adjacency matrix.
#' @export
#'
#' @examples
#' 
part_cor_thres <- function(mat, thres){
  
  p <- dim(mat)[1] 
  
  
  adja <- diag(p)
  
  
  adja[abs(mat) > thres] <- 1
  
  diag(adja) <- 0
  
  return(adja)
  
}
