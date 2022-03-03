
summary.POFD <- function(POFD, plot.Freq = T){
  comp <- apply(!is.na(POFD$POFDs),1,all)
  cat("Grid Size : ",length(POFD$Grid),"\nSample Size : ", nrow(POFD$True.Functions),
      "\n% of complete curves : ",mean(comp),
      "\nSummary of missing sections :\n\n")
  p <- sapply(1:ncol(POFD$POFDs.True.Functions), function(i) sum(length(which(!is.na(POFD$POFDs.True.Functions[,i])))))
  barplot(p, xlab = "t", ylab = "Freq")
  print(summary(rowSums(is.na(POFD$POFDs[!comp,]))/length(POFD$Grid)))
}
