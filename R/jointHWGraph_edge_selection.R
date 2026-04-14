

jointHWGraph_edge_selection <- function(permutations, results, expected_number_of_connections=NULL, FDR= TRUE, target_FDR= 0,group_based_threshold=F){
  
  if(is.null(expected_number_of_connections)){
    expected_number_of_connections <- results$p
  }
  
  taus <- permutations$taus
  
  n_con <- permutations$connection_list
  
  if(group_based_threshold==T){
    n_con_true <- matrix(0, nrow = results$n_groups, ncol = length(taus))
  }
  else{
    n_con_true <- c()
  }
  
  for(i in 1:length(taus)){
    
    thres <- taus[i]
    n_con_time_point  <- c()
    
    for(k in 1:results$n_groups){
      adjacency_matrix <- part_cor_thres(cov2cor(results$omega_list[[k]]), thres = thres)
      
      diag(adjacency_matrix) <- 0
      adjacency_matrix[is.na(adjacency_matrix)] <- 0
      adjacency_matrix[is.nan(adjacency_matrix)] <- 0
      adjacency_matrix[is.null(adjacency_matrix)] <- 0
      
      n_con_time_point[k] <-  (sum(adjacency_matrix)/2) 
      
    }
    
    
    if(group_based_threshold==T){
      n_con_true[,i] <- n_con_time_point
    }
    else{
      n_con_true[i] <- mean(n_con_time_point)
    }
    
  }
  
  FP <- n_con
  
  TP <- n_con_true - n_con
  
  FN <- expected_number_of_connections - TP
  
  F1s <- 2*TP/(2*TP + FP + FN) 
  
  F1s[is.nan(F1s)] <- 0
  
  FDRs <- FP/(FP+TP)
  
  FDRs[is.nan(FDRs)] <- 0
  
  adjacency_matrices <- list()
  
  
  if(FDR == FALSE){
    for(i in 1:results$n_groups){
      
      if(group_based_threshold==T){
        thres <-  taus[F1s[i,]==max(F1s[i,])][1]
      }
      else{
        thres <-  taus[F1s==max(F1s)][1]
      }
      
      adjacency_matrices[[i]]  <- part_cor_thres(cov2cor(results$omega_list[[i]]), thres = thres)
      
    }
    
    if(group_based_threshold==T){
      tau <- c()
      for(i in 1:results$n_groups){
        tau[i] <-  taus[F1s[i,]==max(F1s[i,])][1]
      }
    }
    else{
      tau <-  taus[F1s==max(F1s)][1]
    }
    
  }
  else{
    for(i in 1:results$n_groups){
      
      if(group_based_threshold==T){
        thres <-  taus[FDRs[i,]<=target_FDR][1]
      }
      else{
        thres <-  taus[FDRs<=target_FDR][1]
      }
      
      
      
      adjacency_matrices[[i]]  <- part_cor_thres(cov2cor(results$omega_list[[i]]), thres = thres)
      
    }
    if(group_based_threshold==T){
      tau <- c()
      for(i in 1:results$n_groups){
        tau[i] <-  taus[FDRs[i,]<=target_FDR][1]
      }
    }
    else{
      tau <-  taus[FDRs<=target_FDR][1]
    }
  }
  
  
  
  return(list(adjacency_matrices=adjacency_matrices, tau = tau ) )
}
  
  