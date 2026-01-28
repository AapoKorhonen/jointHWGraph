
threhold_list <- function(taus,map_list,n_timepoints,rep_times){
  
  
  n_con_tau_p <- c()
  
  for(i in 1:length(taus)){
    
    thres <- taus[i]
    
    n_con_tim_p  <- c()
    for(k in 1:n_timepoints){
      
      
      n_con_tim_p_reps <- c()
      for(ii in 1:rep_times){
        adjacency_matrix <- part_cor_thres(cov2cor(map_list[[ii]][[k]]), thres = thres)
        
        diag(adjacency_matrix) <- 0
        
        adjacency_matrix[is.na(adjacency_matrix)] <- 0
        adjacency_matrix[is.nan(adjacency_matrix)] <- 0
        adjacency_matrix[is.null(adjacency_matrix)] <- 0
        n_con_tim_p_reps[ii] <- (sum(adjacency_matrix)/2) 
      }
      
      n_con_tim_p[k] <- mean(n_con_tim_p_reps)
      
    }
    
    
    n_con_tau_p[i] <- mean(n_con_tim_p)
    
  }
  
  return(n_con_tau_p)
}