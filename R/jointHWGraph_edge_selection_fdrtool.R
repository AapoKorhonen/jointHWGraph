
jointHWGraph_edge_selection_fdrtool <- function(jointHWGraph_results, target_FDR = NULL){
  
  p <- jointHWGraph_results$p
  
  
  adjacency_matrices <- list()
  for(i in 1:jointHWGraph_results$time_points){
    
    partial_correlations <- cov2cor(jointHWGraph_results$omega[[i]])
    
    values <- partial_correlations[lower.tri(partial_correlations)]
    z_values <- values
    fdr_tul <- fdrtool::fdrtool(z_values, statistic = "correlation")
    
    
    if(is.null(target_FDR)){
      list11 <- fdr_tul$qval <= fdr_tul$param[1]
    }
    else{
      list11 <- fdr_tul$qval <= target_FDR
    }
    
    list <- rep(0,length(z_values))
    list[list11] <- 1
    
    ad1 <- matrix(0,ncol =  p, nrow =  p)
    ad1[lower.tri(ad1)] <- list
    ad1 <-ad1 +t(ad1) 
    adjacency_matrices[[i]] <- ad1
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
