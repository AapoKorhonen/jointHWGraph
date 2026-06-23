

select_empirical_target_matrix <- function(data, number_or_groups,n, memory_save=F, scale_data=T){
  
  comp_data <- data[[1]]
  
  for(i in 2:number_or_groups){
    comp_data <- rbind(comp_data,data[[i]])
  }
  #if(scale_data) comp_data <- scale(comp_data)
  if(memory_save){ gc()}
  comp_data <- na.omit(comp_data)
  if(memory_save){ gc()}
  est <- cvCovEst::linearShrinkLWEst(comp_data)
  if(memory_save){ gc()}
  est <-chol2inv(chol(est))
  if(memory_save){ gc()}
  B <- diag(diag(est))
  
  if(memory_save){ gc()}
  return(B)
}
