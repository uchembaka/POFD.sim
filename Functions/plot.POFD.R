library(tidyverse)
library(plotly)
library(RColorBrewer)

#Data : a POFD data object
#no.obs : number of observations to include in plot
#plots: type of plots to produce. Combination of any of the following 
  #"noisy.fragments", "true.fragments","obs.grid", "covariance grid", "all.fragmented", "all.true.functions", "true.mean.function", 
  #and "true covariance surface"
# cov.plot.type : "3D" or "2D"
#index: provide specific indexes to plot

#plots = c("noisy.fragments", "true.fragments","obs.grid", "covariance.grid", "all.fragmented", "all.true.functions", "true.mean.function", "true.covariance.surface")

plot.POFD <- function(pofd.obj, nos.obs = 5, plots = c("noisy.fragments", "true.fragments"), cov.plot.type = "3D", index = NA){
  
  cols <- brewer.pal(8, "Set1")
  par(mfrow = c(1,1)); par(mar=c(5, 4, 4, 2) + 0.5)
  
  ### Checks ###
  
  n <- pofd.obj$POFD.args$n
  if(!is.na(index)[1]) {p <- index; nos.obs <- length(p)}
  else p <- sample(n, nos.obs)
  
  if(length(plots) < 1) plots = c("noisy.fragments", "true.fragments")
  
  is.irreg <- pofd.obj$POFD.args$irreg.domain
  if(is.irreg) x <- pofd.obj$Grid$Reg.Grid
  else x <- pofd.obj$Grid
  
  if(pofd.obj$POFD.args$classify) pofd.obj$True.Functions <- pofd.obj$True.Functions$True.Functions


  ########## Plots #############
  y.lims <- c(min(pofd.obj$Dense.Functions), max(pofd.obj$Dense.Functions))
  y.lims2 <- c(min(pofd.obj$True.Functions), max(pofd.obj$True.Functions))
  
  if("true.mean.function" %in% plots){
    if(pofd.obj$POFD.args$classify){
      plot(x, pofd.obj$True.Means$mu.C1,
           type = "l", xlab = "t",
           main = "Mean Functions", lwd=2,
           ylab = expression(paste(mu,"(t)")),
           ylim = y.lims2, cex = 1.5, cex.lab = 1.5)
      lines(x, pofd.obj$True.Means$mu.C2, lwd=2, col = 2)
      legend("bottomright", legend = c(expression(mu[1]), expression(mu[2])), lwd = 2, col = 1:2,
             inset=c(0,1), xpd=TRUE, horiz=TRUE)
    }else{
      plot(x,pofd.obj$True.Mean, type = "l", xlab = "t", main = "Mean Functions",
           ylab = expression(paste(mu,"(t)")), cex = 1.5)
    }
  }
  
  if("true.covariance.surface" %in% plots){
    if(cov.plot.type == "3D"){
      fig <- plot_ly(x =~ x, y =~ x, z=~pofd.obj$True.Covariance,showscale = FALSE)%>%add_surface()%>%add_surface()%>%
        layout(title =  "True Covariance",scene = list(
          xaxis = list(title = "s",titlefont = list(size = 15), tickfont = list(size = 10)),
          yaxis = list(title = "t", titlefont = list(size = 15), tickfont = list(size = 10)),
          zaxis = list(title = TeX("$\\sigma(s,t)$"),titlefont = list(size = 15), tickfont = list(size = 10)),
          margin = list(l = 5, r = 5, b=0, t=0)))%>% config(mathjax = "cdn")
      print(fig)
      print("Check Viewer pane for true Covariance surface")
    }else{
      image(pofd.obj$True.Covariance, cex = 1.5, cex.lab = 1.5, xlab = "t", ylab = "s",
            main = expression(paste(sigma,"(s,t)")))
      contour(pofd.obj$True.Covariance, add = TRUE)
    }
  }
  
  if("all.true.functions" %in% plots){
    matplot(x ,t(pofd.obj$True.Functions),
            type = "l", ylab = "x(t)", xlab = "t", ylim = y.lims2, lty=1,
            main = "True Functions", cex = 1.5, cex.lab = 1.5)
  }
  
  if("all.fragmented" %in% plots){
    matplot(x ,t(pofd.obj$Dense.Functions),
            type = "l", ylab = "y(t)", xlab = "t", ylim = y.lims, lty=1,
            main = "Densely Observed Functions", cex = 1.5, cex.lab = 1.5)
  }
  
  if(!pofd.obj$POFD.args$irreg.domain & ("covariance.grid" %in% plots)){
    tmp <- reshape2::melt(is.na(cov(pofd.obj$POFDs.True.Functions, use="pairwise.complete.obs"))*1)
    tmp$value <- factor(tmp$value)
    plot4 <- ggplot(subset(tmp, value ==0), aes(x = Var2, y = Var1))+
      geom_point()+theme_bw()+labs(x="t", y = "s", title = paste0(pofd.obj$Type, " Covariance")) +
      theme(text = element_text(size = 15, face="bold"))
    print(plot4)
  }
  
  if(!pofd.obj$POFD.args$irreg.domain & ("obs.grid" %in% plots)){
    tmp <- reshape2::melt(is.na(pofd.obj$POFDs.True.Functions)*1)
    tmp$value <- factor(tmp$value)
    plot3 <- ggplot(subset(tmp, value ==0), aes(x = Var2, y = Var1))+
      geom_point()+theme_bw()+labs(x="t", y = "x(t)", title = paste0(pofd.obj$Type, " Observation Grid"))+
      theme(text = element_text(size = 15, face="bold"))
    print(plot3)
  }
  
  if("true.fragments" %in% plots){
    p.lims <- c(min(pofd.obj$True.Functions[p,]), max(pofd.obj$True.Functions[p,]))
    matplot(x, t(pofd.obj$True.Functions[p,, drop = FALSE]), type = "l", lty = 2, col = cols[1:nos.obs],
            ylab = "y(t)", xlab = "t", ylim = p.lims,cex = 1.5, cex.lab = 1.5,
            main = paste0(pofd.obj$POFD.args$POFD.type, " Functional Data (True Functions)"))
    if(is.irreg){
      invisible(sapply(1:length(p), function(i) points(pofd.obj$POFDs.True.Functions$Lt[[p[i]]], pofd.obj$POFDs.True.Functions$Lx[[p[i]]], pch = 19, col = cols[i])))
    }else{
      matpoints(x, t(pofd.obj$POFDs.True.Functions[p,, drop = FALSE]), pch=19, col = cols[1:nos.obs])
    }
  }
  
  if("noisy.fragments" %in% plots){
    p.lims <- c(min(pofd.obj$Dense.Functions[p,]), max(pofd.obj$Dense.Functions[p,]))
    matplot(x, t(pofd.obj$Dense.Functions[p,, drop = FALSE]), type = "l", lty = 2, col = cols[1:nos.obs],
            ylab = "y(t)", xlab = "t", ylim = p.lims,cex = 1.5, cex.lab = 1.5,
            main = paste0(pofd.obj$POFD.args$POFD.type, " Functional Data + Noise"))
    if(is.irreg){
      invisible(sapply(1:length(p), function(i) points(pofd.obj$POFDs$Lt[[p[i]]], pofd.obj$POFDs$Ly[[p[i]]], pch = 19, col = cols[i])))
    }else{
      matpoints(x, t(pofd.obj$POFDs[p,, drop = FALSE]), pch=19, col = cols[1:nos.obs])
    }
  }
  
  
}