
jointHWGraph_permutations <- function(result_object, data, n_permutations=50, EM = T, n_cores=NULL, taus= NULL){
  
  
  n_timepoints <- result_object$time_points
  
  if(is.null(n_cores)){
    n_cores <- parallel::detectCores() - 1
  }
  
  if(is.null(taus)){
    taus <- seq(0, 0.5, length.out=1000)
  }
  
  n_clusters <- n_permutations
  cat("Starting permutations", "\n")
  cl <- parallel::makeCluster(min(n_cores, n_clusters), 
                              type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  
  pb <- progress::progress_bar$new(format = "Permutations :percent [:bar] :elapsed | eta: :eta", 
                                   total = n_permutations + 1, width = 80)
  
  progress <- function() pb$tick()
  pb$tick()
  opts <- list(progress = progress)
  
  `%dopar%` <- foreach::`%dopar%`
  
  permutation_results <- foreach::foreach(i = 1:n_clusters, 
                                          .combine = "rbind", .options.snow = opts) %dopar%{
                              

      suf_data <- data
      for(j in 1:result_object$time_points){
        suf_data[[j]] <- HMFGraph::shuffle_matrix(data[[j]], j)
      }
      
      
       S <- list()
       for(j in 1:n_timepoints){
         S[[j]]  <- cov(suf_data[[j]])
       }

      tvHMFGraph_gibbs <- jointHWGraph_GEM_algorithm(iters=result_object$iters,S = S, p=result_object$p,  n=result_object$n,
                                                     nu= result_object$nu_list,delta = result_object$delta_list,
                                                     epsilon1 = result_object$epsilon1, epsilon2 = result_object$epsilon2,
                                                     fixed_B =  result_object$fixed_B, print_t = F, time_points = n_timepoints,
                                                     print_int=result_object$print_int, stop_criterion=result_object$stop_criterion)
      
      
      map_list <- list()
      m_list <- list()
      for(k in 1:n_timepoints){
        
        m_list[[k]] <- tvHMFGraph_gibbs$omega[[k]]
      }
      
      map_list[[1]] <-  m_list
      
      threhold_list(taus,map_list,n_timepoints=result_object$time_points, rep_times  = 1)
  }
  
  parallel::stopCluster(cl)
  
  return(list(taus = taus , connection_list = colMeans(permutation_results)))
}
