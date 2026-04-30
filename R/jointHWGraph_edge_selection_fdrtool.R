
jointHWGraph_edge_selection_fdrtool <- function(jointHWGraph_results, target_FDR = NULL,verbose=T, plot_fdrtool=F,verbose_fdrtool=F, FDR_results = NULL, memory_save=NULL){
  
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
    z_values <- values
    if(memory_save){
      gc()
    }
    fdr_tul <- fdrtool::fdrtool(z_values, statistic = "correlation",plot=plot_fdrtool,verbose=verbose_fdrtool)
    if(memory_save){
      gc()
    }
    
    
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
    if(memory_save){
      gc()
    }
    if(verbose){
      cat("Network ", i, " ready", "\n")
    }
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
