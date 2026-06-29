#' Title
#'
#' @param data_list 
#' @param delta 
#' @param nu_list 
#' @param n 
#' @param p 
#' @param n_groups 
#' @param iters 
#' @param B 
#' @param stop_criterion 
#' @param scale_data 
#' @param print_t 
#' @param print_int 
#' @param memory_save 
#' @param empirical_B 
#' @param expected_number_of_connections 
#' @param verbose 
#' @param plot_fdrtool 
#' @param verbose_fdrtool 
#' @param limit_FN 
#' @param target_FDR 
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
                         , FDR_control =FALSE , expected_number_of_connections= NULL
                         , verbose=T, plot_fdrtool=F,verbose_fdrtool=F
                         , limit_FN = T, target_FDR = 0.05){
  

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