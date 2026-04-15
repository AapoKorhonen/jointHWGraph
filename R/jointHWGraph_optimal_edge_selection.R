
jointHWGraph_optimal_edge_selection <- function(jointHWGraph_results, expected_number_of_connections= NULL, plot=F,verbose=F){
  
  if(is.null(expected_number_of_connections)){
    expected_number_of_connections <- jointHWGraph_results$p
  }
  
  p <- jointHWGraph_results$p
  
  adjacency_matrices <- list()
  for(i in 1:jointHWGraph_results$n_groups){
    
    partial_correlations <- cov2cor(jointHWGraph_results$omega[[i]])
    
    values <- partial_correlations[lower.tri(partial_correlations)]
    z_values <- values
    fdr_tul <- fdrtool::fdrtool(z_values, statistic = "correlation",plot=plot,verbose=verbose)
    
    
    
    FDRs <- seq(0.00001, 0.99, length.out=1000)
    FP <- c()
    TP <- c()
    
    for(ii in 1:length(FDRs)){
      list11 <- fdr_tul$qval <= FDRs[ii]
      list11[list11 == T] <- 1
      list11[list11 == F] <- 0
      FP[ii] <- floor(FDRs[ii]*sum(list11))
      TP[ii] <- sum(list11) - FP[ii]
      
      if(expected_number_of_connections == 0){
        FP[ii]>0
        break
      }
        
    }
    
    N <- (p^2-p)/2 - expected_number_of_connections
    TN <- N - FP
    FN <- expected_number_of_connections - TP
    
    F1 <- (2 * TP)/(2 * TP + FP + FN)
    
    F1[is.nan(F1)] <- 0
    F1[is.na(F1)] <- 0
    
    best_FDR <- FDRs[which(F1 == max(F1))]
    
    
    list11 <- fdr_tul$qval <= best_FDR[1]
    
    list <- rep(0,length(z_values))
    list[list11] <- 1
    
    ad1 <- matrix(0,ncol =  p, nrow =  p)
    ad1[lower.tri(ad1)] <- list
    ad1 <-ad1 +t(ad1) 
    adjacency_matrices[[i]] <- ad1
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
