library(fdapace)
library(fda)
library(ReconstPoFD)
library(doParallel)
library(foreach)
# rep : number of replications
# comp : compare methods
# use.parallel : use doparallel (recommeded for large replications)
# comp.type : fit, scores or mean.cov (mean and covarince used for curve estimation, .)
# ncores: number of cores to use. Defaults to n-2
# metric: "RMSE", "MSE" 
# data.type : POFD or dense

POFD.sim.comp <- function(..., rep = 10, comp = c("PACE", "Kniep"), use.parallel = TRUE, ncores = NULL, 
                          metric = "RMSE", comp.type = "fit", seed = 1111, data.type = "POFD", method.kniep = 'Error>=0_AlignNO'){
  
  
  "%!in%" <-Negate("%in%")

  get.data <- function(){
    if (use.parallel){
      set.seed(seed, kind = "L'Ecuyer-CMRG")
      data <- foreach(i = 1:rep) %dopar% {
        POFD.sim(...)
      }
    }else{
      data <- vector("list", length = rep)
      for(i in 1:rep){
        set.seed(seed+i)
        data[[i]] <- POFD.sim(...)
      }
    }
    return(data)
  }
  
  get.PACE <- function(){
    if(use.parallel){
      pace <- foreach(i = 1:rep) %dopar% {
        set.seed(seed, kind = "L'Ecuyer-CMRG")
        if(data.type == "POFD"){
          if(is.irreg) sim.list <- data[[i]]$POFDs
          else sim.list <- data[[i]]$POFD.List #POFD.sim2List(data[[i]]$POFDs)
        }else{
          sim.list <- data[[i]]$Dense.List #POFD.sim2List(data[[i]]$Dense.Functions)
        }
        
        FPCA(sim.list[[1]], sim.list[[2]], 
             optns = list(nRegGrid = data[[i]]$POFD.args$grid.size, methodMuCovEst = ifelse(data.type == "dense","cross-sectional","smooth")))
      }
    }else{
      pace <- vector("list", length = rep)
      for(i in 1:rep){
        set.seed(seed+i)
        if(data.type == "POFD"){
          sim.list <- data[[i]]$POFDs#POFD.sim2List(data[[i]]$POFDs)
        }else{
          sim.list <- data[[i]]$Dense.List#POFD.sim2List(data[[i]]$Dense.Functions)
        }
        pace[[i]] <- FPCA(sim.list[[1]], sim.list[[2]], optns = list(nRegGrid = data[[i]]$POFD.args$grid.size))
      }
    }
    
    return(pace)
  }
  
  get.Kraus <- function(){
    if(use.parallel){
      kraus <- foreach(i = 1:rep) %dopar% {
        set.seed(seed, kind = "L'Ecuyer-CMRG")
        if(data.type == "POFD"){
          mat <- t(data[[i]]$POFDs)
        }else{
          mat <- t(data[[i]]$Dense.Functions)
        }
        reconstructKraus(mat)
      }
    }else{
      kraus <- vector("list", length = rep)
      for(i in 1:rep){
        set.seed(seed+i)
        if(data.type == "POFD"){
          mat <- t(data[[i]]$POFDs)
        }else{
          mat <- t(data[[i]]$Dense.Functions)
        }
        kraus[[i]] <- reconstructKraus(mat)
      }
    }
    
    return(kraus)
  }
  
  
  get.Kniep <- function(){
    if(use.parallel){
      kniep <- foreach(i = 1:rep) %dopar% {
        set.seed(seed, kind = "L'Ecuyer-CMRG")
        if(data.type == "POFD"){
          if(is.irreg) sim.list <- data[[i]]$POFDs
          else sim.list <- data[[i]]$POFD.List#POFD.sim2List(data[[i]]$POFDs)
        }else{
          sim.list <- data[[i]]$Dense.List#POFD.sim2List(data[[i]]$Dense.Functions)
        }
        
        reconstructKneipLiebl(Ly = sim.list[[1]], Lu = sim.list[[2]], method = method.kniep, nRegGrid = data[[i]]$POFD.args$grid.size)
      }
    }else{
      kniep <- vector("list", length = rep)
      for(i in 1:rep){
        set.seed(seed+i)
        if(data.type == "POFD"){
          sim.list <- data[[i]]$POFDs#POFD.sim2List(data[[i]]$POFDs)
        }else{
          sim.list <- data[[i]]$Dense.List#POFD.sim2List(data[[i]]$Dense.Functions)
        }
        kniep[[i]] <- reconstructKneipLiebl(Ly = sim.list[[1]], Lu = sim.list[[2]], method = method.kniep, nRegGrid = data[[i]]$POFD.args$grid.size)
      }
    }
    
    return(kniep)
  }
  
  
  get.fit <- function(){
    c <- length(comp)
    fit.res <- vector("list", c)
    ind = 0;  comp.names <- c()
    if("pace" %in% comp){
      ind <- ind+1; comp.names <- c(comp.names,"PACE")
      fit.res[[ind]] <- lapply(pace, function(x) t(fitted(x)))
    }
    if("kraus" %in% comp){
      ind <- ind+1; comp.names <- c(comp.names,"Kraus")
      fit.res[[ind]] <- lapply(kraus, function(x) matrix(unlist(x[['Y_reconst_list']]), ncol=length(x$Y_reconst_list)))
    } 
    if("kniep" %in% comp){
      ind <- ind+1; comp.names <- c(comp.names,"Kniep")
      fit.res[[ind]] <- lapply(kniep, function(x) matrix(unlist(x[['Y_reconst_list']]), ncol=length(x$Y_reconst_list)))
    }
    
    names(fit.res) <- comp.names
    return(fit.res)
  }
  
  get.score <- function(){}#todo
  get.mean.cov <- function(){
    c <- length(comp)
    mean.list <- cov.list <- vector("list", c)
    ind = 0;  comp.names <- c()
    if("pace" %in% comp) {
      ind <- ind+1; comp.names <- c(comp.names,"PACE")
      mean.list[[ind]] <- lapply(pace, function(x) x$mu)
      cov.list[[ind]] <- lapply(pace, function(x) x$fittedCov)
    }
    if("kraus" %in% comp){
      ind <- ind+1; comp.names <- c(comp.names,"Kraus")
      mean.list[[ind]] <- lapply(data, function(x) colMeans(x$POFDs))
      cov.list[[ind]] <- lapply(data, function(x) cov(x$POFDs, use = "pairwise.complete.obs"))
    }
    
    if("kniep" %in% comp){
      ind <- ind+1; comp.names <- c(comp.names,"Kniep")
      
      meancov.list <- lapply(1:rep, function(m){
        KniepMeanCov(type = data.type, POFD.obj = data[[m]], irreg=is.irreg)
      })
      
      mean.list[[ind]] <- lapply(meancov.list, function(mc) mc[[1]])
      cov.list[[ind]] <- lapply(meancov.list, function(mc) mc[[2]])
    }
    
    names(mean.list) <- comp.names; names(cov.list) <- comp.names
    return(list("mean estimate" = mean.list, "cov estimate" = cov.list))
  }
  
  
  
  compare <- function(){
    if("fit" %in% comp.type){
      comp.fit <- comp.fitVar <- matrix(NA, nrow = rep, ncol = length(comp))
      fits <- get.fit()
      
      for(i in 1:rep) {
        for(j in 1:length(fits)){
          if(metric == "RMSE"){
            comp.fit[i,j] <- mean( sapply(1:data[[i]]$POFD.args$n, function (k) Metrics::rmse(data[[i]]$True.Functions[k,], fits[[j]][[i]][,k])) )
            comp.fitVar[i,j] <- var( sapply(1:data[[i]]$POFD.args$n, function (k) Metrics::rmse(data[[i]]$True.Functions[k,], fits[[j]][[i]][,k])) )
          }else if(metric == "MSE"){
            comp.fit[i,j] <- mean( sapply(1:data[[i]]$POFD.args$n, function (k) Metrics::mse(data[[i]]$True.Functions[k,], fits[[j]][[i]][,k])) )
            comp.fitVar[i,j] <- var( sapply(1:data[[i]]$POFD.args$n, function (k) Metrics::mse(data[[i]]$True.Functions[k,], fits[[j]][[i]][,k])) )
          }else{
            stop("specify correct metric")
          }
        }
      }
      colnames(comp.fit) <- colnames(comp.fitVar) <- names(fits)
    }
    
    
    if("mean.cov" %in% comp.type){
      comp.mean <- comp.cov <- matrix(NA, nrow = rep, ncol = length(comp))
      means.covs <- get.mean.cov()
      for(i in 1:rep){
        for (j in 1:length(comp)) {
          if(metric == "RMSE"){
            comp.mean[i,j] <- Metrics::rmse(means.covs[[1]][[j]][[i]], data[[i]]$True.Mean)
            comp.cov[i,j] <- Metrics::rmse(means.covs[[2]][[j]][[i]], data[[i]]$True.Covariance)
          }else if (metric == "MSE"){
            comp.mean[i,j] <- Metrics::mse(means.covs[[1]][[j]][[i]], data[[i]]$True.Mean)
            comp.cov[i,j] <- Metrics::mse(means.covs[[2]][[j]][[i]], data[[i]]$True.Covariance)
          }else{
            stop("specify correct metric")
          }
        }
      }
      colnames(comp.mean) <- colnames(comp.cov) <- names(means.covs[[1]])
    }
    
    
    
    if(length(comp.type) == 1 & comp.type == "fit")comp.list <- list("fit" = comp.fit, "fit.var" = comp.fitVar)
    if(length(comp.type) == 1 & comp.type == "mean.cov") comp.list <- list("mean" = comp.mean, "cov" = comp.cov)
    if(all(c("fit", "mean.cov") %in% comp.type)){
      comp.list <- list("fit" = comp.fit, "fit.var" = comp.fitVar, "mean" = comp.mean, "cov" = comp.cov)
    } 
    return(comp.list)
  }
  
  
  
  ### Checks
  if(use.parallel & all(c("doParallel", "foreach") %!in% (.packages())) == "TRUE") stop("Packages doParallel and foreach not loaded")
  comp <- tolower(comp)
  
  tmp.data <- POFD.sim(...)
  if("kraus" %in% comp & (tmp.data$POFD.args$include.full == F | tmp.data$POFD.args$irreg.domain)){
    warning("data simulations not guaranteed to have complete observation. Removing Kraus from comparison")
    comp <- comp[which(comp != "kniep")]
  }
  
  is.irreg <- tmp.data$POFD.args$irreg.domain
  
  if(length(comp) < 2) stop("Specify more than 1 methods on comp")
  # set up parallel computing 
  if(use.parallel){
    
    if(detectCores() > 3){
      ncores <- ifelse(!is.null(ncores), ncores, detectCores()-1)
      if(ncores > 1){
        registerDoParallel(ncores)
        data <- get.data()
        if("kniep" %in% comp) kniep <- get.Kniep()
        if("pace" %in% comp) pace <- get.PACE()
        if("kraus" %in% comp) kraus <- get.Kraus()
        comp.type = tolower(comp.type)
        comp.res <- compare()
        
      }else{
        use.parallel = FALSE
        print("ncores less than 3 Performing simulation serially")
        data <- get.data()
        if("kniep" %in% comp) kniep <- get.Kniep()
        if("pace" %in% comp) pace <- get.PACE()
        if("kraus" %in% comp) kraus <- get.Kraus()
        comp.type = tolower(comp.type)
        comp.res <- compare()
        
      }
    }
  }else{
    use.parallel = FALSE
    data <- get.data()
    if("kniep" %in% comp) kniep <- get.Kniep()
    if("pace" %in% comp) pace <- get.PACE()
    if("kraus" %in% comp) kraus <- get.Kraus()
    comp.type = tolower(comp.type)
    comp.res <- compare()
  }
  
  out <- list("comp.res"=comp.res, "type" = comp.type, "methods" = paste(comp, collapse = "-"), "rep" = rep, "metric"= metric,
              "POFD.sim.args" = list("n"=data[[1]]$POFD.args$n, "grid.size" = length(data[[1]]$Grid), "error" = data[[1]]$POFD.args$err.sd))
  class(out) <- c("POFD.comp", "list")
  return(out)
  
  
}





irreg2mat <- function(POFD.list, bins){
  
  n <- length(POFD.list[[1]]); bin.size <- length(bins)
  tmp.Ly <- lapply(1:n, function(i) {
    tmp.vec <- rep(NA, bin.size)
    tmp.obs <- which(bins %in% POFD.list$Lt[[i]])
    tmp.vec[tmp.obs] <- POFD.list$Ly[[i]]
    tmp.vec
  })
  res.mat <- matrix(unlist(tmp.Ly), ncol = bin.size, byrow = TRUE)
  
  return(res.mat)
}



KniepMeanCov <- function(type = "POFD", POFD.obj, irreg){
  
  n <- POFD.obj[["POFD.args"]][["n"]]
  
  if(irreg){
    
    #### mean
    
    reg.grid <- POFD.obj[["Grid"]][["Reg.Grid"]]
    if(type == "POFD"){
      argvals <- POFD.obj[["Grid"]][["All.obs.points"]]
      Y <- irreg2mat(POFD.obj$POFDs, argvals)
    }else{
      argvals <- POFD.obj[["Grid"]][["Reg.Grid"]]
      Y <- POFD.obj[["Dense.Functions"]]
    }
    D <- ncol(Y)
    d.vec <- rep(argvals, each = n)
    id <- rep(1:n, rep(D, n))
    gam0 <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = 10))
    mean.full <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    eval.points <- which(argvals %in% reg.grid)
    mean.est <- mean.full[eval.points]
    
    Y.tilde <- Y - matrix(mean.full , n, D, byrow = TRUE)

  }else{
    argvals <- POFD.obj[["Grid"]]
    if(type == "POFD"){
      Y <- POFD.obj[["POFDs"]]
    }else{
      Y <- POFD.obj[["Dense.Functions"]]
    }
    D <- ncol(Y)
    d.vec <- rep(argvals, each = n)
    id <- rep(1:n, rep(D, n))
    gam0 <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = 10))
    mean.est <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    
    Y.tilde <- Y - matrix(mean.est , n, D, byrow = TRUE)
  }
  
  
  ######## Covariance
  
  pev = 0.99
  cov.sum = cov.count = cov.mean = matrix(0,D,D)
  
  for(i in 1:n){
    obs.points = which(!is.na(Y[i,]))
    cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] <- cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0   <- diag(G.0)
  diag(G.0) <- NA
  
  row.vec <- rep(argvals, each = D)
  col.vec <- rep(argvals, D)
  npc.0 <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = 10),weights = as.vector(cov.count)), 
                                    newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  
  # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  w          <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt      <- diag(sqrt(w))
  Winvsqrt   <- diag(1/(sqrt(w)))
  V          <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues    <- replace(evalues, which(evalues <= 0), 0)
  npc        <- length(evalues[evalues>0])
  npc        <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
  efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  # Estimated covariance function
  cov        <- efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  if(irreg) cov <- cov[eval.points, eval.points]
  return(list("mean" = mean.est, "cov" = cov))
}



quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method, 
                trapezoidal = {D <- length(argvals); 1/2 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)], argvals[D] - argvals[D - 1])}, 
                midpoint    = c(0, diff(argvals)), 
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule")
  )
  ##
  return(ret)
}














































































