#' jointHWGraph, EM algorithm
#'
#' @param iters 
#' @param S 
#' @param data_list 
#' @param p 
#' @param n 
#' @param delta 
#' @param nu 
#' @param B 
#' @param print_t 
#' @param n_groups 
#' @param stop_criterion 
#' @param print_int 
#' @param memory_save 
#'
#' @return
#' @export
#'
#' @examples
jointHWGraph_EM_algorithm <- function(iters, S, data_list, p, n, 
                                     delta, nu, B, print_t = TRUE,
                                     n_groups = 1, stop_criterion = 10^(-5),print_int =100,
                                     memory_save = F) {
  
  
  updated_data <- data_list
  
  # Initializing parameters
  
  Phi_0 <- diag(p)
  Omega10 <- diag(p)
  norms <- c()
  L1 <- diag(p)
  shape <- 1
  rate <- 1
  kk <- n_groups 
  B_i <- B
  current_iter_omega <- lapply(seq_len(n_groups), function(i) diag(1, p, p))
  
  missing_groups <- c()
  missing_vals <- list()
  missing <- F
  
  if(memory_save){
    gc()
  }
  # Checking for missing values 
  
  
  for(i in 1:n_groups){
    
    missing_vals[[i]] <- which(is.na(rowSums(data_list[[i]])))
    
    if(length(missing_vals[[i]]) >= 1 ){
      missing <- T
      missing_groups <- append(missing_groups, i)
    }
  }
  
  for (i in seq_len(iters)) {
    
    
    # E-step: Updating missing values
    
    if (missing) {
      for (mi in missing_groups) {
        S[[mi]] <- diag(p)*0
        for (mv in missing_vals[[mi]] ) {
          missing_variables <- which(is.na(data_list[[mi]][mv,
          ]))
          omega_missing <- solve(current_iter_omega[[mi]][missing_variables,
                                                          missing_variables])
          updated_data[[mi]][mv, missing_variables] <- t(-1 *
                                                           omega_missing %*% current_iter_omega[[mi]][missing_variables,
                                                                                                      -missing_variables] %*% data_list[[mi]][mv,
                                                                                                                                              -missing_variables])
          S[[mi]][missing_variables,missing_variables] <- S[[mi]][missing_variables,missing_variables]  + omega_missing
        }
        
        S[[mi]] <- (S[[mi]] +  t(updated_data[[mi]])%*%updated_data[[mi]] )/ n[mi]
        
        if (memory_save) {
          gc()
        }
      }
    }
    
    # E-step: Updating Phi matrix
    
    L1 <- 0
    deg <- 0
    
    for(lk in 1:n_groups){
      L1 <- L1 + (nu[lk] + p - 1)*current_iter_omega[[lk]]
      deg <- deg + nu[lk] + p - 1
    }
    
    L1 <- L1 + (delta + p - 1)*B_i
    
    if(memory_save){
      gc()
    }
    
    L1 <- chol(L1)
    if(memory_save){
      gc()
    }
    
    L1 <- chol2inv(L1)
    if(memory_save){
      gc()
    }
    
    deg <- deg + delta + p - 1
    
    
    Phi_0 <- (deg)*L1
    
    
    L1 <- NULL
    
    if(memory_save){
      gc()
    }
    
    # M-step: Updating presicion matrices
    
    for (k in 1:(n_groups)) {
      
      
      L1 <- (nu[k] + p - 1)*Phi_0 + n[k] * S[[k]]
      
      if(memory_save){
        gc()
      }
      
      L1 <- chol( L1 )
      
      if(memory_save){
        gc()
      }
      
      L1 <- chol2inv( L1)
      
      if(memory_save){
        gc()
      }
      
      Omega10 <- (n[k] + nu[k] +p-1-p-1)*L1
      diff <- current_iter_omega[[k]] - Omega10
      norms[k] <-  sqrt(sum((current_iter_omega[[k]] - Omega10)^2))
      current_iter_omega[[k]] <- Omega10
      
      if(memory_save){
        gc()
      }
      
    }
    Omega10 <- NULL
    L1 <- NULL
    
    
    if(memory_save){
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
