library(RColorBrewer)

plot.POFD.comp <- function(comp.obj){
  
  cols = brewer.pal(4, "Set1")
  
  if(comp.obj$type == "fit"){

    m <- ncol(comp.obj$comp.res[["fit"]])
    boxplot(comp.obj$comp.res[["fit"]], pch = 19, col = cols[1:m],main = paste0("Fit (", comp.obj$metric, ")"))
    #boxplot(comp.obj$comp.res[["fit.var"]], pch = 19, col = cols[1:m], main = "Fit Variance")

  }else if(comp.obj$type == "mean.cov"){
    m <- ncol(comp.obj$comp.res[["mean"]])
    boxplot(comp.obj$comp.res[["mean"]], pch = 19, col = cols[1:m],main = paste0("mean (", comp.obj$metric, ")"))
    boxplot(comp.obj$comp.res[["cov"]], pch = 19, col = cols[1:m], main = "covariance")
  }else{
    stop("invalid comparison type")
  }
  
  
  
  # if(length(comp.obj)>1) comp.obj <- list(comp.obj)
  # cols = brewer.pal(4, "Set1")
  # comp.type <- sapply(1:length(comp.obj), function(i) comp.obj[[i]]$type)
  # comp.meth <- sapply(1:length(comp.obj), function(i) comp.obj[[i]]$methods)
  # comp.metric <- sapply(1:length(comp.obj), function(i) comp.obj[[i]]$metric)
  # nos.rep <- sapply(1:length(comp.obj), function(i) comp.obj[[i]]$rep)
  # if(length(comp.obj)>1){
  #   if(length(unique(comp.type)) != 1) stop("comp.obj list contains POFD.comp of different type")
  #   if(length(unique(comp.meth)) != 1) stop("comp.obj list methods are not consistent")
  #   if(length(unique(comp.metric)) != 1) stop("comp.obj list metrics are not consistent")
  #   if(length(unique(nos.rep)) != 1) warning("comp.obj list number of replications not equal")
  # }
  # 
  # if(unique(comp.type) == "fit"){
  #   if(length(comp.type) == 1){
  #     m <- ncol(comp.obj[[1]]$comp.res[["fit"]])
  #     boxplot(comp.obj[[1]]$comp.res[["fit"]], pch = 19, col = cols[1:m],main = paste0("Fit (", comp.metric, ")"))
  #     boxplot(comp.obj[[1]]$comp.res[["fit.var"]], pch = 19, col = cols[1:m], main = "Fit Variance")
  #   }
  # }
}

