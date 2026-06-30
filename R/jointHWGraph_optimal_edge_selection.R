#' Edge selection for jointHWGraph using an optimal cut-off.
#'
#' @param jointHWGraph_results The output object from the jointHWGraph_EM()-function.
#' @param expected_number_of_connections The number of expected connections in the network. This is used during the edge selection procedure if FDR_control = FALSE. On default, this is set to p, e.i. to be the same as the number of variables.
#' @param verbose Controls if information during the edge selection is printed. On default, verbose = TRUE.
#' @param plot_fdrtool If TRUE, plots the output of fdrtool function. On default, plot_fdrtool = FALSE.
#' @param verbose_fdrtool If TRUE, prints the output of fdrtool function. On default, verbose_fdrtool = FALSE.
#' @param memory_save If TRUE, gc()-function is run after every memory intensive operation. On default, this is set to FALSE. Only recommended to use with large networks p > 10,000.
#' @param limit_FN If TRUE, then the approximated number of FALSE NEGATIVES in the network is limited to be at minimum 0. IF FALSE, FN can become a negative value. On default, limit_FN = TRUE. 
#'
#' @return
#' @export
#'
#' @examples
jointHWGraph_optimal_edge_selection <- function(jointHWGraph_results, expected_number_of_connections= NULL,
                                                verbose=T, plot_fdrtool=F,verbose_fdrtool=F, 
                                                memory_save=NULL, limit_FN = T){
  
  if(verbose){
    cat("Selecting edges for the networks using fdrtool \n")
  }
  
  if(is.null(expected_number_of_connections)){
    expected_number_of_connections <- jointHWGraph_results$p
  }
  
  p <- jointHWGraph_results$p
  
  if(is.null(memory_save)){
    memory_save <- jointHWGraph_results$memory_save 
  }
  
  if(verbose){
    cat("In total ", jointHWGraph_results$n_groups, " networks to be constructured \n")
  }
  
  adjacency_matrices <- list()
  for(i in 1:jointHWGraph_results$n_groups){
    
    if(verbose){
      cat("Constructing the network:", i, "\n")
    }
    
    if(memory_save){
      gc()
    }
    partial_correlations <- cov2cor(jointHWGraph_results$omega[[i]])
    
    values <- partial_correlations[lower.tri(partial_correlations)]
    
    if(memory_save){
      gc()
    }
    fdr_tul <- fdrtool::fdrtool(values, statistic = "correlation",plot=plot_fdrtool,verbose=verbose_fdrtool)
    if(memory_save){
      gc()
    }
    
    q_vals_sort <- sort(fdr_tul$qval)
    
    est_num_con <- floor(q_vals_sort*c(1:length(q_vals_sort)))
    
    FP <- est_num_con
    TP <- c(1:length(q_vals_sort)) - FP
    FDR <- est_num_con/(TP+est_num_con)
    
    N <- (p^2-p)/2 - expected_number_of_connections
    TN <- N - FP
    FN <- expected_number_of_connections - TP
    if(limit_FN){
      FN <- pmax(FN, 0)
    }
    F1 <- (2 * TP)/(2 * TP + FP + FN)

    list11 <- fdr_tul$qval <= sort(fdr_tul$qval)[which(F1==max(F1))[length(which(F1==max(F1)))]]
    
    if(max(TP) == 0){
      list11 <- fdr_tul$qval < sort(fdr_tul$qval)[1]
    }
    
    if(memory_save){
      gc()
    }
    list <- rep(0,length(values))
    list[list11] <- 1
    
    ad1 <- matrix(0,ncol =  p, nrow =  p)
    ad1[lower.tri(ad1)] <- list
    ad1 <-ad1 +t(ad1) 
    adjacency_matrices[[i]] <- ad1
    if(memory_save){
      gc()
    }
    if(verbose){
      cat("Network ", i, " ready", "\n")
    }
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
