library(kableExtra)

summary.POFD.comp <- function(comp.obj){
  cat("Number of replications: ", comp.obj$rep,
      "\nNumber of observations: ", comp.obj$POFD.sim.args$n,
      "\nError sd: ",comp.obj$POFD.sim.args$error,
      "\nMean of Results: \n\n")
  if(comp.obj$type == "fit"){
    df <- rbind((colMeans(comp.obj$comp.res[["fit"]])),
                (colMeans(comp.obj$comp.res[["fit.var"]])))
    rownames(df) <- names(comp.obj$comp.res)
    print(df)
  }else if(comp.obj$type == "mean.cov"){
    df <- rbind((colMeans(comp.obj$comp.res[["mean"]])),
                (colMeans(comp.obj$comp.res[["cov"]])))
    rownames(df) <- names(comp.obj$comp.res)
    print(df)
  }else{
    stop("Invalid comparison type")
  }
  
}
