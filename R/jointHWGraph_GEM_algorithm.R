
jointHWGraph_GEM_algorithm <- function(iters, S, p, n, 
                                     delta, nu,
                                     epsilon1, epsilon2,
                                     fixed_B = FALSE, print_t = TRUE,
                                     time_points = 1, burn_in = 500, stop_criterion = 10^(-5),print_int =100) {
  
  
  
  probs <- numeric(0)
  nut <- list()
  deltat <- list()
  Omega <- list()
  Phi <- list()
  
  Phi_ini <- lapply(seq_len(time_points), function(i) diag(1, p, p))
  Omega_ini <- lapply(seq_len(time_points), function(i) diag(1, p, p))
  
  Phi_0 <- diag(p)
  
  last_iter_Phi <- Phi_ini
  last_iter_omega <- Omega_ini
  
  Omega[[1]] <- Omega_ini
  Phi[[1]] <- Phi_ini
  
  shape <- 1
  rate <- 1
  kk <- time_points 
  
  current_iter_nu <- numeric(0)
  current_iter_delta <- numeric(0)
  current_iter_omega <- Omega_ini
  current_iter_Phi   <- Phi_ini
  current_iter_Phi_B   <- Phi_ini
  
  for (i in seq_len(iters)) {
    
    
    current_iter_nu <- numeric(0)
    current_iter_delta <- numeric(0)
    current_iter_lambda <- numeric(0)
    
    B_i <- diag(1, p, p)
    if (!fixed_B) {
      for (ii in 1:p) {
        shape <- (delta[1] + p - 1) * 0.5 + epsilon1
        rate <-  (delta[1] + p - 1)*(Phi_0[ii, ii] * 0.5) + epsilon2
        B_i[ii, ii] <- (shape-1)/rate
      }
    }
    
    L1 <- (nu[1] + p - 1)*current_iter_omega[[1]]
    deg <- nu[1] + p - 1
    for(lk in 2:time_points){
      L1 <- L1 + (nu[lk] + p - 1)*current_iter_omega[[lk]]
      deg <- deg + nu[lk] + p - 1
    }
    L1 <- chol2inv( chol(L1 + (delta[1] + p - 1)*B_i))
    deg <- deg + delta[1] + p - 1
    
    if (fixed_B == T) {
      Phi_0 <- (deg)*L1
    }
    else{
      Phi_0 <- (deg-p-1)*L1
    }
    
    
    if (time_points > 2) {
      for (k in 1:(time_points)) {
        
        
        
        L2 <- chol2inv( chol( (nu[k] + p - 1)*Phi_0 + n[k] * S[[k]] ))
        Omega10 <- (n[k] + nu[k] +p-1-p-1)*L2
        
        current_iter_omega[[k]] <- Omega10
        
        Omega10_prev <- Omega10
      }
      
    }
    
    
    current_iter_nu[kk] <- nu[kk]
    current_iter_delta[kk] <- delta[kk]
    
    
    if(i > 1){
      
      norms <- c()
      for(iv in 1:time_points){
        
        norms[iv] <- norm(Omega[[iv]] - current_iter_omega[[iv]],type="F")
        
      }
      
      
      for(iv in 1:time_points){
        
        #Phi[[iv]] <- current_iter_Phi[[iv]]
        Omega[[iv]] <- current_iter_omega[[iv]]
        
        
      }
      
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
    for(iv in 1:time_points){
      
      #Phi[[iv]] <- current_iter_Phi[[iv]]
      
      Omega[[iv]] <- current_iter_omega[[iv]]
      
    }
    
  }
  
  return(list(phi = Phi_0,
              omega = Omega))
}
