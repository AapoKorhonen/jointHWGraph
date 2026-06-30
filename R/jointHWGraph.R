#' Bayesian joint Gaussian graphical model using fast EM-algorithm and joint hierarchical Wishart prior.
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
#' @param FDR_control If TRUE, the FDR control procedure is used for the edge selection. On default, this is set to FALSE.
#' @param target_FDR If FDR_control = TRUE, this controls the target FDR value. On default, this is set to 0.05.
#' @param expected_number_of_connections The number of expected connections in the network. This is used during the edge selection procedure if FDR_control = FALSE. On default, this is set to p, e.i. to be the same as the number of variables.
#' @param limit_FN If TRUE, then the approximated number of FALSE NEGATIVES in the network is limited to be at minimum 0. IF FALSE, FN can become a negative value. On default, limit_FN = TRUE. 
#' @param verbose Controls if information during the edge selection is printed. On default, verbose = TRUE.
#' @param plot_fdrtool If TRUE, plots the output of fdrtool function. On default, plot_fdrtool = FALSE.
#' @param verbose_fdrtool If TRUE, prints the output of fdrtool function. On default, verbose_fdrtool = FALSE.
#'
#'
#' @return
#' @export
#'
#' @examples
#' 
jointHWGraph <- function(data_list, delta = NULL, nu_list = NULL,  n=NULL
                         , p=NULL,n_groups=NULL,iters = 10000, B= NULL 
                         , stop_criterion = 10^(-5),scale_data= T, print_t = T
                         , print_int = 100, memory_save=F, empirical_B= F
                         , FDR_control =FALSE, target_FDR = 0.05 , expected_number_of_connections= NULL
                         , verbose=T, plot_fdrtool=F,verbose_fdrtool=F
                         , limit_FN = T){
  

  jointHWGraph_result <- jointHWGraph_EM(data_list = data_list, delta = delta, nu_list = nu_list,  n=n
                                       , p=p,n_groups=n_groups,iters = iters, B= B
                                       , stop_criterion =stop_criterion,scale_data= scale_data, print_t = print_t
                                       , print_int = print_int, memory_save=memory_save) 
  
  if(FDR_control){
    jointHWGraph_edge_selection_results <- jointHWGraph_edge_selection_FDR_control(jointHWGraph_result, target_FDR = target_FDR,verbose=verbose
                                            , plot_fdrtool=plot_fdrtool,verbose_fdrtool=verbose_fdrtool, memory_save=memory_save)
  }
  else{
    jointHWGraph_edge_selection_results <- jointHWGraph_optimal_edge_selection(jointHWGraph_result, expected_number_of_connections= expected_number_of_connections,
                                           verbose=verbose, plot_fdrtool=plot_fdrtool,verbose_fdrtool=verbose_fdrtool, 
                                           memory_save=memory_save, limit_FN = limit_FN)
  }
  
  return(list(omega_list=jointHWGraph_result$omega_list, adjacency_matrices = jointHWGraph_edge_selection_results$adjacency_matrices ))
}