"%!in%" <-Negate("%in%")



POFD.sim <- function(n = 50, grid.size = 100, POFD.type = "fragmented", miss.seg = "S", frag.size = "S", include.full = T,
                     single.frag = F, equal.points = F, nos.points = 5, single.equal.frag = F, range.sparse.obs = c(3,10), 
                     equal.sparse = T,  err.sd = 0.125, base.func = 3, classify = F, mean.fun = "(t-0.5)^2", cov.fun = NA, 
                     full.domain = T, sample.mean = T, sample.cov = T, norm.range = c(1,1,10)){
  if(!is.null(norm.range)){
    if(length(norm.range) != 3 | norm.range[1] %!in% c(0,1)) stop("wrong range.norm specification")
  }else{
    norm.range <- c(2,2,2)
  }
  
  if(base.func == 5){
    if (n < 10) n <- 10
    if (n%%5) n <- n-(n%%5)
  }
  
  
  #for sparse POFD: number of observations in each curve
  no.sparse.obs <- round(runif(ifelse(equal.sparse, 1, n), min(3, range.sparse.obs), max(grid.size/10, range.sparse.obs)))
  #error applied to each obs in each curve (eij)
  err.mat <- matrix(rnorm(n*grid.size,0, err.sd), nrow = n, ncol = grid.size)
  
  range.norm <-function(x,a =0,b =1) ( (x-min(x))/(max(x)-min(x)) )*(b-a)+a # Range normalisation function
  
  sparsePOFD <- function(x.mat,y.mat){#Generate sparse POFD set
    po.y <- y.mat; po.x <- x.mat # y.mat: curves with noise; x.mat: true curves
    for(i in 1:n){
      sparse.obs <- setdiff(1:grid.size, round(runif(ifelse(equal.sparse, no.sparse.obs,no.sparse.obs[i]), 1, grid.size)))
      po.y[i,sparse.obs] <- NA; po.x[i,sparse.obs] <- NA
    }
    return(list("po.y"=po.y, "po.x"=po.x))
  }
  
  fragments <- function(){# Generate missing fragments
    
    
    segs <- round(ifelse(single.frag, 1, runif(1, 1, min(5, grid.size/10)))) # generate segments of curve 
    vec<- rep(T, grid.size)
    tj <- 0
    for(s in 1:segs){
      if(equal.points){
        points <- nos.points
      }else{
        points <- switch (frag.size,
                          "S" = round(runif(1, min(3, ceiling(0.3*(grid.size/segs))), ceiling(0.3*(grid.size/segs)) )),
                          "M" = round(runif(1, min(3, ceiling(0.3*(grid.size/segs))), ceiling(0.9*(grid.size/segs)) )),
                          "L" = round(runif(1, ceiling(0.5*(grid.size/segs)), ceiling(1*(grid.size/segs))))
        )
        #points <- round(runif(1, ifelse(frag.size == "S", 3, 9), ceiling((ifelse(frag.size == "S", 0.3, 0.7)*(grid.size/max(3,segs) )))))
      }
      if(s == 1){
        grid.seg <- round(grid.size/segs)
        tj <- ifelse(points>=grid.seg, 1, sample.int((grid.seg-points), 1)) #start of obs point
        end <- tj+points# end of obs point
        
      }else{
        grid.seg <- round(grid.size/segs)
        next.seg <- ifelse(s==segs, grid.size, grid.seg * s)
        prev.seg <- next.seg - grid.seg
        prev.end <- end
        ctrl.tj <- (prev.seg - end)
        tj <- ifelse(points>=grid.seg, next.seg-points, round( runif(1, prev.end, next.seg-points) ))
        end <- tj+points
      }
      vec[tj:end] <- F
    }
    return(vec)
  }
  
  
  complete.domain <- function(x.mat, y.mat, po.x, po.y, min.freq = 1){# ensure sample set covers all time points
    ind = T
    while (ind == T) {
      rnd.curves <- sample.int(nrow(x.mat), round(0.1*nrow(x.mat)))
      nos.used <- 0
      for (i in 1:length(rnd.curves)) {
        col.count <- sapply(1:ncol(x.mat), function(i) sum(length(which(!is.na(po.x[,i])))))
        empty.col <- which (col.count < min.freq)
        if(length(empty.col) == 0){
          ind = F
          break
        } 
        points <- round(runif(1, min(3,0.3*(grid.size/3)), 0.3*(grid.size/3)))
        po.x[rnd.curves[i],] <- NA; po.y[rnd.curves[i],] <- NA
        if( (points + empty.col[1]) > grid.size ){
          po.x[rnd.curves[i], (grid.size-points):grid.size] <- x.mat[rnd.curves[i], (grid.size-points):grid.size]
          po.y[rnd.curves[i], (grid.size-points):grid.size] <- y.mat[rnd.curves[i], (grid.size-points):grid.size]
        }else{
          po.x[rnd.curves[i], empty.col[1]:(points + empty.col[1])] <- x.mat[rnd.curves[i], empty.col[1]:(points + empty.col[1])]
          po.x[rnd.curves[i], empty.col[1]:(points + empty.col[1])] <- x.mat[rnd.curves[i], empty.col[1]:(points + empty.col[1])]
        }
      }
    }
    
    return(list("po.y"=po.y, "po.x"=po.x))
  }
  
  
  fragPOFD <- function(x.mat, y.mat){ #apply generated fragments to full curve
    # y.mat: curves with noise; x.mat: true curves
    po.x <- x.mat; po.y = y.mat
    nos.edge.curves <- runif(1, ifelse(n < 20,1,4), round(ifelse(single.frag, 0.15, 0.10)*n)) #number of curves with frags at beginning and end
    edge.curves.ind <- sample.int(n,nos.edge.curves)# randomly select curves from sample set
    side.count <- 0
    
    for (i in 1:n) {
      if(single.equal.frag){# to be removed
        start = sample(1:(grid.size - (grid.size*0.3)),1)
        end = start+grid.size*0.3
        po.x[i,start:end] <- NA; po.y[i,start:end] <- NA
      }else{
        if(frag.size == "S" & i %in% edge.curves.ind){
          side.count <- side.count+1
          fp <- fragments(edge.ctr = T, ind = ifelse(side.count%%2 == 0, 0,1))
        }else{
          fp <- fragments()
        }
        po.x[i,fp] <- NA
        po.y[i,fp] <- NA
      }
    }
    if(include.full){# include full curve if True
      full.curves <-  sample(1:n, max(1,round((sample(1:20,1)/100)*n)))
      po.x[full.curves,] <- x.mat[full.curves,]
      po.y[full.curves,] <- y.mat[full.curves,]
    }
    
    if(frag.size == "S" & full.domain){
      return(complete.domain(x.mat, y.mat, po.x, po.y))
    }
    
    return(list("po.y"=po.y, "po.x"=po.x))
  }
  
  
  Kraus <- function(){
    grid <- seq(0,1, length.out = grid.size)
    K = 100
    x <- sapply(1:n, function(i){
      ksi <- sqrt(exp(-(((1:K)-1)^2)/5))*rnorm(K)
      eta <- sqrt(exp(-((1:K)^2)/5))*rnorm(K)
      vk = 3^(-(2*(1:K)-1)); wk = 3^(-(2*(1:K)))
      colSums(t(sapply(1:K, function(k) (sqrt(2)*sqrt(vk[k])*ksi[k]*(cos(2*k*t(pi*grid))/sqrt(5))) + (sqrt(2)*sqrt(wk[k])*eta[k]*(sin(2*k*t(pi*grid))/sqrt(5))) )))
    })
    #y <- t(sapply(1:n, function(i) t(x)[i,] + err.mat[i,])); 
    x <-t(x)
    if(norm.range[1] == 1) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    y <- x + err.mat
    
    if(POFD.type == "sparse"){
      pofd <- sparsePOFD(x.mat=x, y.mat=y)
    }else{
      pofd <- fragPOFD(x.mat=x, y.mat=y)
    }
    mu <- colMeans(x)
    return(list("Grid" = grid, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "True.Covariance" = cov(x),
                "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x, "Type" = POFD.type))
  }
  
  
  Alois <- function(mu.f = 1){
    grid <- seq(0,1, length.out = grid.size)
    if(mu.f == 1){
      mu <- grid^(ifelse(base.func == 2, 2,1))+sin(2*pi*grid)
    }else{
      mu <- grid^(ifelse(base.func == 2, 2,1))+cos(2*pi*grid)
    }
    
    K = 50
    variations <- sapply(1:n, function(i){
      ksi.1 <- 50*sqrt(exp(-(((1:K)-1)^2)/ifelse(base.func == 2,1,5)))*rnorm(K)
      ksi.2 <- 50*sqrt(exp(-((1:K)^2)/ifelse(base.func == 2,1,5)))*rnorm(K)  
      #t((t(ksi.1)%*%(cos((1:K)%*%t(pi*grid))/sqrt(5))) + (t(ksi.2)%*%(sin((1:K)%*%t(pi*grid))/sqrt(5))))
      colMeans(t(sapply(1:K, function(k) (ksi.1[k]*(cos(k*t(pi*grid))/sqrt(5))) + (ksi.2[k]*(sin(k*t(pi*grid))/sqrt(5))) )))
    })
    x <- t(sapply(1:n, function(i) mu + t(variations)[i,] ))
    #y <- t(sapply(1:n, function(i) mu + t(variations)[i,] + err.mat[i,]))
    if(norm.range[1] == 1) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    y <- x + err.mat
    
    if(POFD.type == "sparse"){
      pofd <- sparsePOFD(x.mat=x, y.mat=y)
    }else{
      pofd <- fragPOFD(x.mat=x, y.mat=y)
    }
    
    if(sample.mean) mu <- colMeans(x)
    
    return(list("Grid" = grid, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "True.Covariance" = cov(x),
                "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x, "Type" = POFD.type))
  }
  
  
  Wei <- function(Z = runif(n, -2/3,2/3)){
    grid <- 1:grid.size
    mu <- (sin(grid/15)/((grid-50)^2 / 100 + 1)) + (sqrt((3*grid + 100)*(1))/50)
    #y <- t(sapply(1:n, function(i) (sin(grid/15)/((grid-50)^2 / 100 + 1)) + Z[i]*sqrt((3*grid + 100)*(Z[i]+1))/50 + err.mat[i,]))
    x <-  t(sapply(1:n, function(i) (sin(grid/15)/((grid-50)^2 / 100 + 1)) + Z[i]*sqrt((3*grid + 100)*(Z[i]+1))/50))
    if(norm.range[1] == 1) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    y <- x + err.mat
    if(POFD.type == "sparse"){
      pofd <- sparsePOFD(x.mat=x, y.mat=y)
    }else{
      pofd <- fragPOFD(x.mat=x, y.mat=y)
    }
    if(sample.mean) mu <- colMeans(x)
    
    return(list("Grid" = grid, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "True.Covariance" = cov(x),
                "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x, "Type" = POFD.type))
  }
  
  GP <- function(){
    kern <- function(s, t) {
      if(cov.fun =="Wei-1"){
        Wei.kern_1(s,t)
      }else if(!is.na(cov.fun)){
        eval(parse(text = cov.fun))
      }else{
        phis <- list(phi.1 = function(t) 1 + t*0, phi.2 = function(t) (2*t - 1)*sqrt(3),
                     phi.3 <- function(t) (6*t^2 - 6*t + 1)*sqrt(5), phi.4 <- function(t) (20 * t^3 - 30*t^2 +12*t -1)*sqrt(7))
        sum(sapply(1:4, function(i) 0.5^(i-1)*phis[[i]](s)*phis[[i]](t)))
      }
    }# covariance function
    
    mu <- function (t){
      if(is.na(mean.fun)) t*0
      else eval(parse(text = mean.fun))
    } # mean function
    t <- seq(0, 1, length.out = grid.size) # will sample the GP at these points
    Sig <- matrix(nrow = grid.size, ncol = grid.size)
    for (i in 1:grid.size) for (j in 1:grid.size) Sig[i, j] = kern(t[i], t[j])
    x <- MASS::mvrnorm(n, mu(t), Sig)
    if(norm.range[1] == 1) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    y <- x + err.mat
    if(POFD.type == "sparse"){
      pofd <- sparsePOFD(x.mat=x, y.mat=y)
    }else{
      pofd <- fragPOFD(x.mat=x, y.mat=y)
    }
    
    if(sample.mean) mu <- colMeans(x)
    else mu <- mu(t)
    if(sample.cov) cov.x <- cov(x)
    else cov.x <- Sig
    
    return(list("Grid" = t, "True.Mean" = mu, "True.Functions" = x, 
                "Dense.Functions" = y, "True.Covariance" = cov.x,
                "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x, "Type" = POFD.type))
  }
  
  
  composite.functions <- function(){
    t <- seq(0,1, length.out = grid.size)
    x <- matrix(NA, nrow = n, ncol = grid.size)
    z <- runif(n, -2/3, 2/3)
    f1 <- t+sin(2*pi*t)
    f2 <- 7*(t-0.5)^3
    f3 <- 0.5*exp(t)+ (t-0.5)^2
    f4 <- t-cos(2*pi*t)
    f5 <- sin(2*pi*t)+cos(2*pi*t)
    m <- n/5
    for (i in 1:m) {
      j <- (i-1) * 4
      x[j+i,] <- f1*z[j+i]
      x[j+i+1,] <- f2*z[j+i+1]
      x[j+i+2,] <- f3*z[j+i+2]
      x[j+i+3,] <- f4*z[j+i+3]
      x[j+i+4,] <- f5*z[j+i+4]
    }
    if(norm.range[1] == 1) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    y = x+err.mat
    if(POFD.type == "sparse"){
      pofd <- sparsePOFD(x.mat=x, y.mat=y)
    }else{
      pofd <- fragPOFD(x.mat=x, y.mat=y)
    }
    
    mu <- colMeans(x)
    return(list("Grid" = t, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "True.Covariance" = cov(x),
                "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x, "Type" = POFD.type))
    
  }
  
  Wei.kern_1 <- function(s,t){
    phi.sin <- function(k,t) sqrt(2)*sin(2*k*pi*t)
    phi.cos <- function(k,t) sqrt(2)*cos(2*k*pi*t)
    K <- 0
    for (i in 1:4) {
      if(i%%2 == 1){
        K <- K + 0.5^(i-1)*phi.sin(1,s)*phi.sin(1,t)
      }else{
        K <- K + 0.5^(i-1)*phi.cos(2,s)*phi.cos(2,t)
      }
    }
    return(K)
  }
  
  
  makeClass <- function(){# generate class
    if(base.func > 2){
      print("base.func set to 1")
      base.func = 1
    }
    samps <- list(sample(1:n, ceiling(n/2)), sample(1:n, floor(n/2)))
    k <- c(rep(1, length(samps[[1]])),  rep(2, length(samps[[2]])))
    class.1 <- Alois(); class.2 <- Alois(mu.f = 2)
    
    mu1 <- class.1$True.Mean; x1 <- class.1$True.Functions[samps[[1]],]
    mu2 <- class.2$True.Mean; x2 <- class.2$True.Functions[samps[[2]],]
    y <- rbind(class.1$Dense.Functions[samps[[1]],], class.2$Dense.Functions[samps[[2]],])
    x <- rbind(x1,x2)
    pofd.y <- rbind(class.1$POFDs[samps[[1]],], class.2$POFDs[samps[[2]],])
    pofd.x <- rbind(class.1$POFDs.True.Functions[samps[[1]],], class.2$POFDs.True.Functions[samps[[2]],])
    return(list("Grid" = class.1$Grid, "True.Functions" = x, "True.Mean.1" = mu1, "True.Functions.1" = x1, 
                "True.Mean.2" = mu2, "True.Functions.2" = x2, "Dense.Functions" = y, "True.Covariance" = cov(x), 
                "POFDs" = pofd.y, "POFDs.True.Functions" = pofd.x, "classes" = k, "Type" = POFD.type))
  }
  
  
  if(classify){
    sim <- makeClass()
    class(sim) <- c("POFD", "list")
    return(sim)
  }else{
    if(base.func %in% c(1,2)){
      sim <- Alois()
      class(sim) <- c("POFD", "list")
      return(sim)
    }else if(base.func == 3){
      sim <- Kraus()
      class(sim) <- c("POFD", "list")
      return(sim)
    }else if(base.func == 4){
      sim <- Wei()
      class(sim) <- c("POFD", "list")
      return(sim)
    }else if (base.func == 5){
      sim <- composite.functions()
      class(sim) <- c("POFD", "list")
      return(sim)
    }else{
      sim <- GP()
      class(sim) <- c("POFD", "list")
      return(sim)
    }
  }
  
  
}#POFD.sim

#Data : a POFD data object
#no.obs : number of observations to include in plot
#with.class : for POFD with class
#no.plots: number of plots to produce. {2: POFD plots with and without noise
#3: }
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


summary.POFD <- function(POFD, plot.Freq = T){
  comp <- apply(!is.na(POFD$POFDs),1,all)
  cat("Grid Size : ",length(POFD$Grid),"\nSample Size : ", nrow(POFD$True.Functions),
      "\n% of complete curves : ",mean(comp),
      "\nSummary of missing sections :\n\n")
  p <- sapply(1:ncol(POFD$POFDs.True.Functions), function(i) sum(length(which(!is.na(POFD$POFDs.True.Functions[,i])))))
  barplot(p, xlab = "t", ylab = "Freq")
  print(summary(rowSums(is.na(POFD$POFDs[!comp,]))/length(POFD$Grid)))
}

POFD.sim2List <- function(FD.mat){
  grid <- seq(0,1, len=ncol(FD.mat)) 
  Ly <- split(t(FD.mat), rep(1:ncol(t(FD.mat)), each = nrow(t(FD.mat))))
  Lt <- lapply(Ly, function(i) grid[!is.na(i)])
  Ly <- lapply(Ly, function(i) i[!is.na(i)])
  return(list(Ly, Lt))
}



myFunction <- function(data){
  if(any(class(data) == "POFD")) data <- data$POFDs
  else if(any(class(data) == "matrix")) data <- data
  else stop("Data not of class matrix or POFD")
  newdata <- data
  for (i in 1:nrow(data)){
    tmp.mean.vec <- colMeans(data[-i,], na.rm = T)
    newdata <- rbind(newdata, tmp.mean.vec)
  }
  
  return(newdata)
}



Kniep2POFD.sim <- function(simuldata){
  n <- ncol(simuldata$U_true_mat)
  grid.size <- nrow(simuldata$U_true_mat)
  m <- nrow(simuldata$U_mat)
  grid <- simuldata$U_true_mat[,1]
  mu <- simuldata$mean_true_vec
  x <- t(simuldata$Y_true_mat)
  err.mat <- matrix(rnorm(n*grid.size,0, 0.0125), nrow = n, ncol = grid.size)
  y <- x+err.mat
  obs.grid <- sapply(1:n, function (i) {
    sapply(1:m, function(j){
      birk::which.closest(round(grid,3), round(simuldata$U_mat[j,i],3))
    })
  })
  po.x <- po.y <- matrix(NA, nrow = n, ncol = grid.size)
  
  for(i in 1:n){
    po.y[i, obs.grid[,i]] <- simuldata$Y_mat[,i]
    po.x[i, obs.grid[,i]] <- x[i,obs.grid[,i]]
    y[i, obs.grid[,i]] <- simuldata$Y_mat[,i]
  } 
  
  sim <- list("Grid" = grid, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "True.Covariance" = cov(x),
              "POFDs" = po.y, "POFDs.True.Functions" = po.x, "Type" = "Fragmented")
  class(sim) <- c("POFD", "list")
  return(sim)
}


## TO be included in POFD.sim
Mercer.kern.decomp <- function(s, t, k=1:4, lambda = "0.5^(k-1)", phi = list(phi.sin = "sqrt(2)*sin(2*k*t)", 
                                                                                     phi.cos = "sqrt(2)*cos(2*k*t)"), alternate.k = T, repeat.phi = T){
  
  if(length(phi)%%2 != 0) stop("Incorrect number of basis functions specified")
  if(!repeat.phi & length(k) != length(phi)) stop("Number of basis functions not equal lenght of K and repeat.phi set to F")
  
  if(length(k) == length(phi)){
    phi.list <- vector(mode = "list", length = length(k))
    for(i in 1:length(k)){
      phi.list[[i]] <- function(t) eval(parse(text = phi[[i]]))
    }
    return(sum(sapply(k, function(i) eval(parse(text = lambda))*phi.list[[i]](s)*phi.list[[i]](t))))
  }
  
  if(length(k) != length(phi) & repeat.phi){
    if(length(phi) > 2) stop("Number of basis function must be 2, i.e sin(kt) and cos(kt)")
    phi.list <- vector(mode = "list", length = 2)
    for(i in 1:length(phi)){
      phi.list[[i]] <- function(t) eval(parse(text = phi[[i]]))
    }
    if(alternate.k){
      phi1 <- sapply(which(k%%2 == 1), function(i) eval(parse(text = lambda))*phi.list[[1]](s)*phi.list[[1]](t))
      phi2 <- sapply(which(k%%2 == 0), function(i) eval(parse(text = lambda))*phi.list[[2]](s)*phi.list[[2]](t))
      return(sum(phi1,phi2))
    }else{
      return(sum(sapply(k, function(i) eval(parse(text = lambda))*phi.list[[i]](s)*phi.list[[i]](t))))
    }
  }
}