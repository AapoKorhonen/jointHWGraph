
jointHWGraph_GEM_algorithm <- function(iters, S, data_list, p, n, 
                                     delta, nu,
                                     epsilon1, epsilon2,
                                     fixed_B = FALSE, print_t = TRUE,
                                     time_points = 1, burn_in = 500, stop_criterion = 10^(-5),print_int =100) {
  
  
  updated_data <- data_list
  
  
  Phi_0 <- diag(p)
  Omega10 <- diag(p)
  norms <- c()
  L1 <- diag(p)
  shape <- 1
  rate <- 1
  kk <- time_points 
  B_i <- diag(1, p, p)
  current_iter_omega <- lapply(seq_len(time_points), function(i) diag(1, p, p))
  
  missing_groups <- c()
  missing_vals <- list()
  missing <- F
  
  for(i in 1:time_points){
    
    missing_vals[[i]] <- which(is.na(rowSums(data_list[[i]])))
    
    if(length(missing_vals[[i]]) >= 1 ){
      missing <- T
      missing_groups <- append(missing_groups, i)
    }
  }
  
  for (i in seq_len(iters)) {
    
    
    if (!fixed_B) {
      for (ii in 1:p) {
        shape <- (delta[1] + p - 1) * 0.5 + epsilon1
        rate <-  (delta[1] + p - 1)*(Phi_0[ii, ii] * 0.5) + epsilon2
        B_i[ii, ii] <- (shape-1)/rate
      }
    }
    
    if(missing){
      for(mi in missing_groups){
        
        for(mv in missing_vals[[mi]]){
          
          missing_variables <- which(is.na(data_list[[mi]][mv,]))
          
          updated_data[[mi]][mv,missing_variables] <- t(-1*solve(current_iter_omega[[mi]][ missing_variables, missing_variables])
                                                 %*%current_iter_omega[[mi]][ missing_variables, -missing_variables]
                                                 %*%data_list[[mi]][mv,-missing_variables])
          
        }
        
        S[[mi]] <- cov(updated_data[[mi]])
      }
      
    }
    
    L1 <- (nu[1] + p - 1)*current_iter_omega[[1]]
    deg <- nu[1] + p - 1
    for(lk in 2:time_points){
      L1 <- L1 + (nu[lk] + p - 1)*current_iter_omega[[lk]]
      deg <- deg + nu[lk] + p - 1
    }
    L1 <- L1 + (delta[1] + p - 1)*B_i
    L1 <- chol(L1)
    L1 <- chol2inv(L1)
    deg <- deg + delta[1] + p - 1
    
    if (fixed_B == T) {
      Phi_0[] <- (deg)*L1
    }
    else{
      Phi_0[] <- (deg-p-1)*L1
    }
    
    L1 <- NULL
    
    if(p>=1000){
      gc()
    }
    
    for (k in 1:(time_points)) {
      
      L1 <- (nu[k] + p - 1)*Phi_0 + n[k] * S[[k]]
      L1 <- chol( L1 )
      L1 <- chol2inv( L1)
      Omega10 <- (n[k] + nu[k] +p-1-p-1)*L1
      diff <- current_iter_omega[[k]] - Omega10
      norms[k] <-  sqrt(sum((current_iter_omega[[k]] - Omega10)^2))
      #norms[k] <- norm(current_iter_omega[[k]] - Omega10,type="F")
      current_iter_omega[[k]] <- Omega10
    }
    Omega10 <- NULL
    L1 <- NULL
    if(p>=10000){
      gc()
    }
    
    if(i > 1){
      
      if(print_t==T){
        if(i %% print_int == 0){
          print(i)
          print(mean(norms))
        }
      }
      
      if(mean(norms) < stop_criterion){
        if(print_t==T){
          print(i)
          print("stop")
        }
        break
      }
    }
    
  }
  
  return(list(phi = Phi_0,
              omega = current_iter_omega))
}
