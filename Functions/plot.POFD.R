
#Data : a POFD data object
#no.obs : number of observations to include in plot
#with.class : for POFD with class
#no.plots: number of plots to produce.

plot.POFD <- function(data, no.obs = 5, with.class = F, no.plots = 2, index = NA){
  #if (length(dev.list()!=0)) {dev.off()}
  library(tidyverse); library(plotly)
  
  y.lims <- c(min(data$Dense.Functions), max(data$Dense.Functions))
  y.lims2 <- c(min(data$True.Functions), max(data$True.Functions))
  n <- nrow(data$Dense.Functions)
  if(is.na(index)) p <- sample(n, no.obs)
  else p <- index
  no.plots <- ifelse(no.plots < 1, 1, no.plots)
  x <- data$Grid;
  
  if(no.plots >= 8){
    
    if(with.class){
      plot(x,data$True.Mean.1, type = "l", xlab = "t", main = "Mean Functions", lwd=2,
           ylab = expression(paste(mu,"(t)")), ylim = c(min(data$Dense.Functions), max(data$Dense.Functions)))
      lines(x, data$True.Mean.2, lwd=2, col = 2)
      legend("topleft", legend = c("Function 1", "Function 2"), lwd = 2, col = 1:2)
    }else{
      plot(x,data$True.Mean, type = "l", xlab = "t", main = "Mean Functions",
           ylab = expression(paste(mu,"(t)")))
    }
  }
  if(no.plots >= 8){
    fig <- plot_ly(z=~data$True.Covariance)%>%add_surface()
    print(fig)
    print("Check Viewer pane for true Covariance surface")
  }
  if(no.plots >= 6){
    matplot(x ,t(data$Dense.Functions), type = "l", ylab = "y(t)", xlab = "t", ylim = y.lims, lty=1, main = "Densely Observed Functions")
  }
  if(no.plots >= 5){
    matplot(x ,t(data$True.Functions), type = "l", ylab = "x(t)", xlab = "t", ylim = y.lims2, lty=1, main = "True Functions")
  }
  if(no.plots >= 4){
    tmp <- reshape2::melt(is.na(cov(data$POFDs.True.Functions, use="pairwise.complete.obs"))*1)
    tmp$value <- factor(tmp$value)
    plot4 <- ggplot(subset(tmp, value ==0), aes(x = Var2, y = Var1))+geom_point()+theme_bw()+labs(x="t", y = "s", title = paste0(data$Type, " Covariance"))
    print(plot4)
  }
  if(no.plots >= 3){
    tmp <- reshape2::melt(is.na(data$POFDs.True.Functions)*1)
    tmp$value <- factor(tmp$value)
    plot3 <- ggplot(subset(tmp, value ==0), aes(x = Var2, y = Var1))+geom_point()+theme_bw()+labs(x="t", y = "x(t)", title = paste0(data$Type, " Observation Grid"))
    print(plot3)
  }
  if(no.plots >= 2){
    p.lims2 <- c(min(data$True.Functions[p,]), max(data$True.Functions[p,]))
    matplot(x ,t(data$POFDs.True.Functions[p,, drop = FALSE]), pch=19, type = "b", ylab = "x(t)", xlab = "t", ylim = p.lims2, lty=1, main = paste0(data$Type, " Functional Data (True Functions)"))
    matlines(x ,t(data$True.Functions[p,, drop = FALSE]), type = "l", lty=2, col = 1)
  }
  p.lims <- c(min(data$Dense.Functions[p,]), max(data$Dense.Functions[p,]))
  matplot(x ,t(data$POFDs[p,, drop = FALSE]), pch=19, type = "b", ylab = "y(t)", xlab = "t", ylim = p.lims, lty=1, main = paste0(data$Type, " Functional Data + Noise"))
  matlines(x ,t(data$Dense.Functions[p,, drop = FALSE]), type = "l", lty=2, col = 1)
  
}#plot