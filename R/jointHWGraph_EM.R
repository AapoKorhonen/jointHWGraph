#' EM algorithm for jointHWGraph
#'
#' @param data_list a list of data matrices
#' @param delta value for delta parameter. Default value is p*10
#' @param nu_list vector of nu values. If a single value is given, then it is used for all groups. Default value is p*10
#' @param n vector of sample sizes
#' @param p number of variables 
#' @param n_groups number of groups
#' @param iters max iterations used for the gem algorithm
#' @param B The target matrix. On default, an identity matrix is used 
#' @param stop_criterion Stopping criterion. Default value is 10^-5
#' @param scale_data If TRUE, then data is scale with scale-function
#' @param print_t If TRUE, then
#' @param only_mean 
#' @param print_int 
#' @param memory_save This option should only be used with extremely large networks p>10,000. If TRUE, then gc-function is used after all matrix operations to save memory
#'
#' @return
#' @export
#'
#' @examples
#' 
jointHWGraph_EM <- function(data_list, delta = NULL, nu_list = NULL,  n=NULL
                            , p=NULL,n_groups=NULL,iters = 10000, B= NULL 
                            , stop_criterion = 10^(-5),scale_data= T, print_t = T
                            , print_int = 100, memory_save=F){
  
  if(is.null(n_groups)){
    n_groups = length(data_list)
  }
  if(is.null(p)){
    p = dim(data_list[[1]])[2]
  }
  
  if(is.null(n)){
    n <- c()
    for(i in 1:n_groups){
      n[i] <- dim(data_list[[i]])[1]
    }
    
  }
  
  if(is.null(delta)){
    delta <- p*10
  }
  
  if(is.null(nu_list)){
    nu_list <- rep(p*10,n_groups)
  }
  else if (length(nu_list)==1){
    
    nu_list <- rep(nu_list,n_groups)
    
  }
  
  S <- list()
  for(i in 1:n_groups){
    if(scale_data){
      data_list[[i]] <- scale(data_list[[i]])
    }
    S[[i]]  <- (t(data_list[[i]])%*%data_list[[i]])/(n[i])
  }
  
  if(is.null(B)){
    B <- diag(p)
  }
  
  if(memory_save){
    gc()
  }
  
  jointHWGraph_result <- jointHWGraph_EM_algorithm(iters=iters,S = S, data_list=data_list, p=p,  n=n, 
                                               nu= nu_list,delta = delta, B = B,
                                               print_t = print_t, n_groups = n_groups,
                                               print_int=print_int, stop=stop_criterion,
                                               memory_save = memory_save) 
  
  
  return(list(omega_list = jointHWGraph_result$omega,phi = jointHWGraph_result$phi, iters=iters,
              p = p, n = n, nu_list = nu_list, delta= delta, n_groups = n_groups,
              print_int=print_int, stop_criterion=stop_criterion, memory_save = memory_save))
  
}
