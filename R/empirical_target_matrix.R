

select_empirical_target_matrix <- function(data, number_or_groups,n, memory_save=F, cut_off = 1){
  
  comp_data <- data[[1]]
  
  for(i in 2:number_or_groups){
    comp_data <- rbind(comp_data,data[[i]])
  }
  
  #t(comp_data)%*%comp_data/(sum(n))
  if(memory_save){ gc()}
  est <- cvCovEst::linearShrinkLWEst(comp_data)
  if(memory_save){ gc()}
  est <-chol2inv(chol(est))
  if(memory_save){ gc()}
  B <- diag(diag(est))
  
  if(max(diag(B))/min(diag(B)) < cut_off ) B <- diag(p)
  if(memory_save){ gc()}
  return(B)
}