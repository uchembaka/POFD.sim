
summary.POFD <- function(POFD.obj, plot.Freq = T){
  
  if(POFD.obj$POFD.args$irreg.domain){
    nos.obs <- unlist(lapply(POFD.obj$POFDs$Ly, function(y) length(y)))
    cat("Grid Size : ",POFD.obj$POFD.args$grid.size,"\nSample Size : ",POFD.obj$POFD.args$n,
        "\nMin number of obs: ", min(nos.obs),
        "\nMax number of obs: ", max(nos.obs),
        "\nMedian number of obs: ", median(nos.obs), "\n")
  }else{
    comp <- apply(!is.na(POFD.obj$POFDs),1,all)
    cat("Grid Size : ",POFD.obj$POFD.args$grid.size,"\nSample Size : ",POFD.obj$POFD.args$n ,
        "\n% of complete curves : ",mean(comp),
        "\nSummary of missing sections :\n\n")
    p <- sapply(1:ncol(POFD.obj$POFDs.True.Functions), function(i) sum(length(which(!is.na(POFD.obj$POFDs.True.Functions[,i])))))
    barplot(p, xlab = "t", ylab = "Freq")
    print(summary(rowSums(is.na(POFD.obj$POFDs[!comp,]))/POFD.obj$POFD.args$grid.size))
  }

}
