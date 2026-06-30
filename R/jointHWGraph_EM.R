#' EM algorithm for jointHWGraph
#'
#' @param data_list A list of data matrices n x p. 
#' @param delta Hyperparameter value. On default value of 1 is used. Controls how much  shrinkage is used on off-diagonal elements. 
#' @param nu_list A list, or single value for hyperparameter values. If single value is given, then that is used for all groups. On default nu = p*10 for all groups. Controls how much strongly groups are shrank to each other. 
#' @param n Number of sample. List or a single value. On default, this is derived from the data_list.
#' @param p Number of variables. On default, this is derived from the data_list.
#' @param n_groups Number of groups. On default, this is derived from the data_list.
#' @param iters Number of max iters for the EM algorithm.
#' @param B Hyperparameter value. This should be a positive-definite matrix. On default, identity matrix is used.
#' @param stop_criterion Convergence criterion. On default, 10^-5 is used.
#' @param scale_data If TRUE, the data is scaled to have mean zero and variance 1. On default, this is TRUE.
#' @param print_t If TRUE, then convergence information is printed during the EM algorithm
#' @param print_int Defines how often convergence information is printed during the EM algorithm. On default 100 is used, e.i. information is printed every 100 iterations.
#' @param memory_save If TRUE, gc()-function is run after every memory intensive operation. On default, this is set to FALSE. Only recommended to use with large networks p > 10,000.
#' @param empirical_B If TRUE, an empirical B matrix is calculated with linearShrinkLWEst()-function. On default, this is set to FALSE.
#'
#' @return
#' @export
#'
#' @examples
#' 
jointHWGraph_EM <- function(data_list, delta = NULL, nu_list = NULL,  n=NULL
                            , p=NULL,n_groups=NULL,iters = 10000, B= NULL 
                            , stop_criterion = 10^(-5),scale_data= T, print_t = T
                            , print_int = 100, memory_save=F, empirical_B= F){
  
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
    delta <- 1
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
  
  if(empirical_B){
    B <- select_empirical_target_matrix(data_list, number_or_groups=n_groups,n=n, memory_save=memory_save,scale_data=scale_data)
  }
  
  if(memory_save){
    gc()
  }
  
  jointHWGraph_result <- jointHWGraph_EM_algorithm(iters=iters,S = S, data_list=data_list, p=p,  n=n, 
                                               nu= nu_list,delta = delta, B = B,
                                               print_t = print_t, n_groups = n_groups,
                                               print_int=print_int, stop=stop_criterion,
                                               memory_save = memory_save) 
  
  
  return(list(omega_list = jointHWGraph_result$omega, iters=iters,
              p = p, n = n, nu_list = nu_list, delta= delta, n_groups = n_groups,
              print_int=print_int, stop_criterion=stop_criterion, memory_save = memory_save))
  
}
