
jointHWGraph_difference_networks_fdrtool <- function(jointHWGraph_results, target_FDR = NULL, plot=F,verbose=F){
  
  p <- jointHWGraph_results$p
  
  
  adjacency_matrices <- list()
  
  for(i in 1:(jointHWGraph_results$n_groups-1)){
    
    adjacency_matrices_per <- list()
    
    partial_correlations1 <- cov2cor(jointHWGraph_results$omega[[i]])
    
    partial_correlations1_Z <- GeneNet::z.transform(partial_correlations1)
    
    for(ii in (i+1):jointHWGraph_results$n_groups ){
      
      
      partial_correlations2 <- cov2cor(jointHWGraph_results$omega[[ii]])
      
      partial_correlations2_Z <- GeneNet::z.transform(partial_correlations2)
      
      differences <- partial_correlations1_Z - partial_correlations2_Z
      
      values <- differences[lower.tri(differences)]
      
      fdr_tul <- fdrtool::fdrtool(values, statistic = "correlation",plot=plot,verbose=verbose)
      
      
      if(is.null(target_FDR)){
        list11 <- fdr_tul$qval <= fdr_tul$param[1]
      }
      else{
        list11 <- fdr_tul$qval <= target_FDR
      }
      
      list <- rep(0,length(values))
      list[list11] <- 1
      
      ad1 <- matrix(0,ncol =  p, nrow =  p)
      ad1[lower.tri(ad1)] <- list
      ad1 <-ad1 +t(ad1) 
      adjacency_matrices_per[[ii - i  ]] <- ad1
    }
    adjacency_matrices[[i]] <- adjacency_matrices_per
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
