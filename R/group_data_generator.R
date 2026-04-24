#' Samples multiple datasets with different random network partial correlation structures. Similarity of the networks is controllable.
#'
#' @param n The number of samples for different groups. A list or a single value. If a single value, then all groups have the same number of samples. 
#' @param p The number of variables.
#' @param d The number of connection in the network. Default value is p.
#' @param similarity The similarity between groups' networks. Default value is 0.50
#' @param number_of_groups The number of groups. Default value is 2.
#'
#' @return 
#' @export
#'
#' @example path.R
group_data_generator <- function(n, p, d, similarity = 0.50, number_of_groups=2,positives=0.5){
  
  
  
  
  if(d==0){
    d = p
  }
  
  if(length(n) == 1){
    n <- rep(n,number_of_groups)
  }
  
  
  # Firstly the graph structure is generated
  
  presicion_matrix <- replicate(number_of_groups,diag(p))
                           
  common <- floor(similarity*d)
  
  differences <- d - common
  
  com_network <- igraph::sample_gnm(n = p,m= (common+differences*number_of_groups) ,directed = FALSE)
  
  connections_S <- rep(0, (common+differences*number_of_groups))
  
  for(i in 1:(common+differences*number_of_groups)){
    if (positives < runif(1)){
      connections_S[i] <- 1
    }
    else{
      connections_S[i] <- -1
      
      
    }
  }
  
  
  edges <- igraph::as_adjacency_matrix(com_network)[lower.tri(igraph::as_adjacency_matrix(com_network))]

  edges[edges==1] <- connections_S
  
  
  mat <- diag(p)*0
  
  mat[lower.tri(mat)] <- edges
  
  mat <- mat+ t(mat)
  
  list_differences <- (common+1):(common+differences*number_of_groups)
  
  presicion_list <- list()
  for(i in 1:number_of_groups){
    
    
    removed <- c(1:differences)+(i-1)*differences
    presicion_list[[i]] <- mat*as.matrix(igraph::as_adjacency_matrix(igraph::delete_edges(com_network,list_differences[-1*removed] )))
    
  }
  
  omega_list <- list()
  sigma_list <- list()
  data_list <- list()
  adjacency_list <- list()
  
  for(i in 1:number_of_groups){
    
    adjacency_list[[i]] <- abs(presicion_list[[i]])
    
    omega <- presicion_list[[i]] 
    
    diag(omega) <- 0
    
    diag(omega) = abs(min(eigen(omega,symmetric=T,only.values=T)$values)) + 0.1

    sigma = cov2cor(chol2inv(chol(omega)))
    sigma_list[[i]] <- sigma
    omega = chol2inv(chol(sigma))
    
    omega_list[[i]] <- omega
    # Using Rcpp for sampling. This is much faster than mvrnorm in R if p >> 100.
    
    x = mvrnorm_cpp(n[i], rep(0, p), sigma)
    
    data_list[[i]] <- x

  }
  
  
  return(list(data = data_list, precision_matrix = omega_list, covariance_matrix = sigma_list, adjacency_matrix=adjacency_list))
}

