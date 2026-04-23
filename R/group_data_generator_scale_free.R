#' Generate group data with similar scale-free partial correlation structure
#'
#' @param n 
#' @param p 
#' @param number_of_groups 
#' @param positives 
#' @param similarity 
#'
#' @return
#' @export
#'
#' @examples
group_data_generator_scale_free <- function(n=100, p=100 ,number_of_groups=2,positives=0.5, similarity=0.9){
  
  
  if(length(n) == 1){
    n <- rep(n,number_of_groups)
  }
  
  
  adjacency_matrices <- list()
  
  org_adjacency_matrices <- diag(p)
  diag(org_adjacency_matrices) <- 0
  
  if(runif(1) > positives){# half of the partial correlations are set to be positive and the other half negatives
    org_adjacency_matrices[1,2] <- 1
    org_adjacency_matrices[2,1] <- 1
  }
  else{
    org_adjacency_matrices[1,2] <- -1
    org_adjacency_matrices[2,1] <- -1
  }
  
  for(i in 3:p){ 
    degrees <- colSums(abs(org_adjacency_matrices))
    con <- sample(1:p,1 , prob =  degrees/sum(degrees) )
    
    # half of the partial correlations are set to be positive and the other half negatives
    if(runif(1) > positives){
      org_adjacency_matrices[i,con] <- 1
      org_adjacency_matrices[con,i] <- org_adjacency_matrices[i,con]
    }
    else{
      org_adjacency_matrices[i,con] <- -1
      org_adjacency_matrices[con,i] <- org_adjacency_matrices[i,con]
    }
  }
  
  diag(org_adjacency_matrices) <- 0
  
  part_cor_mats <- list()
  x <- list()
  sigma <- list()
  omega <- list()
  
  unique_connections <- floor((1-similarity)*(p-1))
  
  if(unique_connections>0){
    for(ii in 1:unique_connections){
      # Remove connections from the "original/common" network
      
      csum <- colSums(abs(org_adjacency_matrices) )
      values <- which(csum>0)
      ind <- sample(1:p,1,prob = csum/sum(csum)  )
      ind <- sample(values,1 )
      
      numbers <- which(abs(org_adjacency_matrices[ind,])==1)
      if(length(numbers)==1){
        ind2 <- numbers
      }
      else{
        ind2 <- sample(numbers,1)
      }
      org_adjacency_matrices[ind,ind2] <- 0
      org_adjacency_matrices[ind2,ind] <- 0
    }
  }
  
  diag(org_adjacency_matrices) <- 1 # set diagonals to value 1 
  for (i in 1:(number_of_groups)){
    adjacency_matrices[[i]] <- org_adjacency_matrices
  }  
  
  if(unique_connections>0){
    
    for(ii in 1:(unique_connections)){
      
      
      for(i in 1:number_of_groups){
        
        csum <- colSums(abs(adjacency_matrices[[i]]))
        
        probs <- (csum-1) # degree + 1, because diagonal elements are 1 
        #ind3 <- sample(1:p,1 )
        if(sum(probs)==0){
          ind3 <- 1
        }
        else{
          ind3 <- sample(1:p,1,prob= probs/sum(probs) )
        }
        
        
        arvot2 <- which(org_adjacency_matrices[ind3,] == 0) # check the common and other networks, selects an unique connection
        arvot <- arvot2
        
        if(length(arvot)==1){
          ind4 <- arvot
        }
        else{
          #ind4 <- sample(arvot,1,prob = (probs/sum(probs))[arvot])
          ind4 <- sample(arvot,1)
        }
        
        if(runif(1) > positives){ # half of the partial correlations are set to be positive and the other half negatives
          adjacency_matrices[[i]][ind3,ind4] <- 1
          adjacency_matrices[[i]][ind4,ind3] <- 1
          org_adjacency_matrices[ind3,ind4] <- 1 # include the new connection to org_adjacency for easy checking
          org_adjacency_matrices[ind4,ind3] <- 1 
        }
        else{
          adjacency_matrices[[i]][ind3,ind4] <- -1
          adjacency_matrices[[i]][ind4,ind3] <- -1
          org_adjacency_matrices[ind3,ind4] <- -1 # include the new connection to org_adjacency for easy checking
          org_adjacency_matrices[ind4,ind3] <- -1 
        }
        
        
      }
      
    }
  }
  
  omega_list <- list()
  sigma_list <- list()
  data_list <- list()
  
  for (i in 1:number_of_groups) {
    
    diag(adjacency_matrices[[i]]) = 0
    
    omega[[i]] = adjacency_matrices[[i]]
    
    diag(omega[[i]]) = abs(min(eigen(omega[[i]] , only.values = T, symmetric = T )$values)) + 0.1
    
    sigma_list[[i]] = cov2cor(   chol2inv(  chol(omega[[i]])  ))
    omega_list[[i]] = chol2inv(  chol(sigma_list[[i]]))
    
    x = HMFGraph::mvrnorm_cpp(n[i], rep(0, p), sigma_list[[i]])
    data_list[[i]] <- x
    adjacency_matrices[[i]] <- abs(adjacency_matrices[[i]])
  }
  
  sim = list(data = data_list, covariance_matrix = sigma_list, precision_matrix = omega_list, adjacency_matrix= adjacency_matrices)
  
  return(sim)
}
