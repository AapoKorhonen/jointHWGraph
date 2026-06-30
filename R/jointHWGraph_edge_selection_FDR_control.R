#' Edge selection for jointHWGraph using FDR control.
#'
#' @param jointHWGraph_results The output object from the jointHWGraph_EM()-function.
#' @param target_FDR If FDR_control = TRUE, this controls the target FDR value. On default, this is set to 0.05.
#' @param verbose Controls if information during the edge selection is printed. On default, verbose = TRUE.
#' @param plot_fdrtool If TRUE, plots the output of fdrtool function. On default, plot_fdrtool = FALSE.
#' @param verbose_fdrtool If TRUE, prints the output of fdrtool function. On default, verbose_fdrtool = FALSE.
#' @param memory_save If TRUE, gc()-function is run after every memory intensive operation. On default, this is set to FALSE. Only recommended to use with large networks p > 10,000.
#'
#' @return
#' @export
#'
#' @examples
#' 
jointHWGraph_edge_selection_FDR_control <- function(jointHWGraph_results, target_FDR = 0.05,verbose=T, plot_fdrtool=F,verbose_fdrtool=F, memory_save=NULL){
  
  p <- jointHWGraph_results$p
  
  if(verbose){
    cat("Selecting edges for the networks using fdrtool \n")
  }
  
  if(is.null(memory_save)){
    memory_save <- jointHWGraph_results$memory_save 
  }
  
  adjacency_matrices <- list()
  if(verbose){
    cat("In total ", jointHWGraph_results$n_groups, " networks to be constructured \n")
  }
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
    
    list11 <- fdr_tul$qval <= target_FDR
  
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
