
jointHWGraph_edge_selection_HC <- function(jointHWGraph_results,limit=F, alpha0=0.1,verbose=T, plot_fdrtool=F,verbose_fdrtool=F, FDR_results = NULL, memory_save=NULL){
  
  p <- jointHWGraph_results$p
  
  # Functions from fdrtool package. These are not visible outside the package and thus
  # we had to copy them here in order to have a efficient computation
  F0 = function(x, param) {
    return(fdrtool::pcor0(x, kappa = param))
  }
  get.pval = function(x, param) {
    ax = abs(x)
    return(ifelse(ax == Inf, 0, pmax(.Machine$double.eps, 
                                     2 - 2 * F0(ax, param))))
  }
  
  if(verbose){
    cat("Selecting edges for the networks using fdrtool \n")
  }
  
  if(is.null(memory_save)){
    memory_save <- jointHWGraph_results$memory_save 
  }
  
  adjacency_matrices <- list()
  if(verbose){
    cat("In total ", jointHWGraph_results$n_groups, " networks to be constructured \n")
  }
  for(i in 1:jointHWGraph_results$n_groups){
    if(verbose){
      cat("Constructing the network:", i, "\n")
    }
    if(memory_save){
      gc()
    }
    partial_correlations <- cov2cor(jointHWGraph_results$omega[[i]])
    
    values <- partial_correlations[lower.tri(partial_correlations)]
    
    if(memory_save){
      gc()
    }
    
    x0 =  fdrtool::fndr.cutoff(values, "correlation")
    if(memory_save){
      gc()
    }
    
    cf.out <- fdrtool::censored.fit(x = values, cutoff = x0, statistic = "correlation")
    scale.param <- cf.out[1, 5]
    
    if(memory_save){
      gc()
    }
    
    pval = get.pval(values, scale.param)
    
    if(memory_save){
      gc()
    }
    
    hc_res <- fdrtool::hc.thresh(pval, alpha0 = alpha0,plot = plot_fdrtool)
    
    
    if(memory_save){
      gc()
    }
    list11 <- pval <= hc_res
    
    
    list <- rep(0,length(values))
    list[list11] <- 1
    
    ad1 <- matrix(0,ncol =  p, nrow =  p)
    ad1[lower.tri(ad1)] <- list
    ad1 <-ad1 +t(ad1) 
    adjacency_matrices[[i]] <- ad1
    if(memory_save){
      gc()
    }
    if(verbose){
      cat("Network ", i, " ready", "\n")
    }
  }
  
  
  return(list(adjacency_matrices=adjacency_matrices))
}
