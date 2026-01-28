
jointHWGraph_GEM <- function(data_list, delta_list = NULL, nu_list = NULL,  n=NULL
                            , p=NULL,time_points=NULL,iters = 10000, n_networks = NULL, B= diag(p) 
                            , stop_criterion = 10^(-5)
                            , inter=100,epsilon1 = 0.001
                            , epsilon2 = 0.001,fixed_B = T, print_t = T, w=0.5
                            , fixed_nu = T, fixed_delta = T, sd = 10
                            , burn_in=500, only_mean=F, print_int = 100){
  
  if(is.null(time_points)){
    time_points = length(data_list)
  }
  if(is.null(p)){
    p = dim(data_list[[1]])[2]
    
    
  }
  if(is.null(n)){
    n <- c()
    for(i in 1:time_points){
      n[i] <- dim(data_list[[i]])[1]
    }
    
  }
  
  if(is.null(delta_list)){
    delta_list <- rep(p*10,time_points)
  }
  else if (length(delta_list)==1){
    
    delta_list <- rep(delta_list,time_points)
    
  }
  if(is.null(nu_list)){
    nu_list <- rep(p*10,time_points)
  }
  else if (length(nu_list)==1){
    
    nu_list <- rep(nu_list,time_points)
    
  }
  
  S <- list()
  for(i in 1:time_points){
    S[[i]]  <- cov(data_list[[i]])
  }
  print(p)
  
  tvHMFGraph_result <- jointHWGraph_GEM_algorithm(iters=iters,S = S, p=p,  n=n, 
                                               nu= nu_list,delta = delta_list, epsilon1 = epsilon1, 
                                               epsilon2 = epsilon2,fixed_B =fixed_B, 
                                               print_t = print_t, time_points = time_points,
                                               print_int=print_int, stop=stop_criterion) 
  
  
  return(list(omega_list = tvHMFGraph_result$omega,iters=iters,
              p = p, n = n, nu_list = nu_list, delta_list= delta_list, time_points = time_points,
              print_int=print_int, stop_criterion=stop_criterion, epsilon1 = epsilon1, epsilon2 = epsilon2,
              fixed_B  = fixed_B))
  
}
