## Version 3.0
library(boot)


# n: number of functions in sample
# grid.size: For regular grid simlation. Number of points if the fully observed curves
# POFD.type: Type of the POFD - "fragmented" or "sparse"
# frag.size: Type of fragment ("S" for short, "M" for mixed and "L" for long)
# include.full: Specify weather to include fully observed curves
# single.frag: Specify whether each curve should have oonly one fragment per curve (fragmented curves only)
# equal.points: Whether each curve should have equal number of observed points
# nos.points: Number of points to use if equal.points = TRUE
# sinle.equal.frag: NO LONGER IN USE
# range.sparse.obs: Range of number of observation in the sparese functional data
# irreg.domain: Whether to use irreg domain for fragmented samples
# equal.sparse: Whether each curve should have equal number of points in the sparse setting (TO BE REMOVED)
# err.sd: Standard deviation of zero mean normal error added to each obs
# base.func; A 2 item list with first item specifying the base function used for simulation and second argument for base 
#function specifications. Also accept a single number corresponding the base function. (Default settings will apply).
#Base function 1 & 2: Base function specified in Kneip and Liebl (2020) (No optn.args needed)
#Base function 3: BAse function specified in Kraus (2015) (No optn.args needed)
#Base function 4: Gaussian process. with mean and covariance function specified using mean.fun and cov.fun
#Base function 5: equal composition of functions of the form t+sin(2*pi*t), 7*(t-0.5)^3, 0.5*exp(t)+(t-0.5)^2, t-cos(2*pi*t)
# sin(2*pi*t) + cos(2*pi*t)
#Base function 6: Base functions specified in Delaigle and Hall (2013) (optns.args xt and mu see paper for ref)

# classify: Whether to make create class of function (uses base.fun = 6)
# mean.fun: Mean specification for GP functions should be specified in terms of t. e.g "(t-0.5)^2"
# cov.fun: Covariance specification for GP. should be specified in terms of s and t. e.g "0.005*exp((-abs(s-t)^2/0.05))". One can also
#specify covariance as matern (use args.matern to specify args or NULL for default) and mercer's decomposition use args.mercer to specify form
# Full.domain: Whether the fragmented observations in short case should be guaranteed to cover the entire domain.
# norm.range: normalize the range of the function. The first entry of the vector should be 1 or 0 specifing whether or not to 
#normalize range (1 for TRUE) to the second (min) and third (max) entry.
# args.Mercer: specify decomposition of mercer. args are K: vector of index, lambda: form of eigen values e.g "0.5^(k-1)",
#phi: list of eigen function e.g "list("sqrt(2)*sin(2*k*t)", "sqrt(2)*cos(2*k*t)"), alternate.k, whether the phi's should be alternated based on length(k)
#repeat.phi: repeat the eigen function as function of k. See default for example

# args.matern: specify shape and range args for matern e.g list(l = 1, v = 1)


POFD.sim <- function(n = 50, grid.size = 100, grid.range = c(0,1), POFD.type = "fragmented", miss.seg = "S", frag.size = "S", include.full = TRUE,
                     single.frag = FALSE, equal.points = FALSE, nos.points = 5, single.equal.frag = FALSE, range.sparse.obs = c(3,15), irreg.domain = FALSE,
                     equal.sparse = TRUE,  err.sd = 0.125, base.func = list(func=1, optn.args = NULL), classify = FALSE, mean.fun = "(t-0.5)^2", cov.fun = "0.005*exp((-abs(s-t)^2/0.05))", 
                     full.domain = TRUE, norm.range = c(0,0,10),
                     args.Mercer = list(k = 1:4, lambda = "0.5^(k-1)", phi = list(phi.sin = "sqrt(2)*sin(2*k*t)",phi.cos = "sqrt(2)*cos(2*k*t)"), alternate.k = T, repeat.phi = T),
                     Matern.args = list(l = 1, v = 1)){
  
  "%!in%" <-Negate("%in%")
  
  ############################################# Checks #####################################################
  if(class(base.func) == "numeric"){
    base.func <- list(func = base.func, optn.args = NULL)
  }
  
  if(norm.range[1] == 1){
    change.rng <- TRUE
    if (length(norm.range) != 3) stop("norm.range not correctly specified")
  }else{
    change.rng <- FALSE
  }
  
  if(base.func[[1]] > 6){
    warning("Value specified for base.func[[\"func\"]] is > 6, setting to 1")
    base.func[[1]] = 1
  } 
  
  if(base.func[[1]] == 5){
    if (n < 10) n <- 10
    if (n%%5) n <- n-(n%%5)
  }
  
  if(base.func[[1]] == 6 & is.null(base.func[[2]])){
    print("Using default optional arguments")
    base.func$optn.args <- list(xt = 1, mu = 1)
  }
  
  if(classify & base.func[[1]] != 6){
    base.func[[1]] = 6
    print("base.func$func set to 6")
    if(is.null(base.func[[2]])){
      print("Using default optional arguments")
      base.func[[2]] <- list(xt = 1, mu = 1)
    }
  }
  if(base.func[[1]] == 6 & grid.range[2] < 100) grid.range <- c(0,100)
  
  frag.size <- toupper(frag.size)
  if(frag.size %!in% c("S", "M", "L")) stop("invalid fragment size. Use \"S\" or \"M\" or \"L\" ")
  POFD.type <- tolower(POFD.type)
  
  if(base.func[[1]] != 4){
    cov.fun = NULL
    mean.fun = NULL
    args.Mercer = NULL
    Matern.args = NULL
  }else if (is.character(cov.fun)){
    if(tolower(cov.fun) == "mercer") Matern.args = NULL
    if(tolower(cov.fun) == "matern") args.Mercer = NULL
    if(tolower(cov.fun) != "matern" & tolower(cov.fun) != "mercer"){
      args.Mercer = NULL
      Matern.args = NULL
    }
  }else{
    args.Mercer = NULL
    Matern.args = NULL
  }
  
  
  call.args <- as.list(environment())
  
  #for sparse POFD: number of observations in each curve
  no.sparse.obs <- round(runif(ifelse(equal.sparse, 1, n), min(3, range.sparse.obs), max(grid.size/10, range.sparse.obs)))
  
  
  ######################################### Inner functions #########################################################
  
  
  range.norm <-function(x,a =0,b =1) ( (x-min(x))/(max(x)-min(x)) )*(b-a)+a # Range normalisation function
  
  list2mat <- function(dataList){
    Ls <- lapply(dataList, function(dl) dl[1:grid.size])
    mat <- matrix(unlist(Ls), ncol = grid.size, byrow = TRUE)
    return(mat)
  }
  
  mat2List <- function(FD.mat){
    grid <- seq(grid.range[1], grid.range[2], len=ncol(FD.mat)) 
    Ly <- split(t(FD.mat), rep(1:ncol(t(FD.mat)), each = nrow(t(FD.mat))))
    Lt <- lapply(Ly, function(i) grid[!is.na(i)])
    Ly <- lapply(Ly, function(i) i[!is.na(i)])
    return(list("Ly" = Ly,"Lt" = Lt))
  }
  
  sparsePOFD <- function(x.mat,y.mat){#Generate sparse POFD set
    po.y <- y.mat; po.x <- x.mat # y.mat: curves with noise; x.mat: true curves
    for(i in 1:n){
      sparse.obs <- setdiff(1:grid.size, round(runif(ifelse(equal.sparse, no.sparse.obs,no.sparse.obs[i]), 1, grid.size)))
      po.y[i,sparse.obs] <- NA; po.x[i,sparse.obs] <- NA
    }
    return(list("po.y"=po.y, "po.x"=po.x))
  }
  
  fragments <- function(){# Generate missing fragments
    
    segs <- round(ifelse(single.frag, 1, runif(1, 1, min(3, grid.size/10)))) # generate segments of curve 
    vec<- rep(TRUE, grid.size)
    tj <- 0
    seg.total <- segs*(floor(grid.size/segs))#total available points based on segments partitions
    off <- grid.size - seg.total
    for(s in 1:segs){
      if(equal.points){
        points <- nos.points
      }else{
        points <- switch (frag.size,
                          "S" = round(runif(1, min(3, ceiling(0.2*(grid.size/segs))), ceiling(0.2*(grid.size/segs)) )),
                          "M" = round(runif(1, min(4, ceiling(0.3*(grid.size/segs))), ifelse(s==segs, floor(1*(grid.size/segs))+off, floor(1*(grid.size/segs))) )),
                          "L" = round(runif(1, ceiling(0.6*(grid.size/segs)), ifelse(s==segs, floor(1*(grid.size/segs))+off, floor(1*(grid.size/segs)))  ))
        )
      }
      if(s == 1){
        grid.seg <- floor(grid.size/segs)
        tj <- ifelse(points>=grid.seg, 1, sample.int((grid.seg-points), 1)) #start of obs point
        end <- tj+points-1# end of obs point
        
      }else{
        grid.seg <- floor(grid.size/segs)
        next.seg <- ifelse(s==segs, grid.size, grid.seg * s)
        prev.seg <- next.seg - grid.seg
        prev.end <- end
        ctrl.tj <- (prev.seg - end)
        tj <- ifelse(points>=grid.seg, next.seg-points, floor( runif(1, prev.end, next.seg-points) ))
        end <- tj+points-1
      }
      vec[tj:end] <- F
    }
    return(vec)
  }
  
  
  
  irregFragment <- function(){
    
    t.reg <- seq(grid.range[1],grid.range[2], length.out = grid.size)
    t.frag.list <- vector("list", n)
    for (i in 1: n){
      segs <- ifelse(single.frag, 1, sample.int(3,1))
      t.frag <- round(t.reg,4)
      seg.total <- segs*(floor(grid.size/segs))#total available points based on segments partitions
      off <- grid.size - seg.total
      for (s in 1:segs) {
        if(equal.points){
          points <- nos.points
        }else{
          points <- switch (frag.size,
                            "S" = round(runif(1, min(3, ceiling(0.3*(grid.size/segs))), ceiling(0.3*(grid.size/segs)) )),
                            "M" = round(runif(1, min(4, ceiling(0.3*(grid.size/segs))), ifelse(s==segs, floor(1*(grid.size/segs))+off, floor(1*(grid.size/segs))) )),
                            "L" = round(runif(1, ceiling(0.6*(grid.size/segs)), ifelse(s==segs, floor(1*(grid.size/segs))+off, floor(1*(grid.size/segs)))  ))
          )
        }
        if(s == 1){
          if(segs == 1){
            ai <- sample(c(0,0.3,0.6), 1)
            bi <- ifelse(ai == 0.6, 1, ai+0.3)
          }else{
            ai <- 0; bi <- 0.3
          }
        }else if (s ==2){
          if(segs == 2){
            ai <- sample(c(0.3,0.6), 1)
            bi <- ifelse(ai == 0.6, 1, ai+0.3)
          }
        }else{
          ai <- 0.6; bi <- 1
        }
        
        t.frag <- c(t.frag, round(sort(unique(runif(points, ai*grid.range[2], bi*grid.range[2]))),4))
      }
      
      t.frag.list[[i]] <- sort(unique(t.frag))
    }
    
    
    return(t.frag.list)
    
  }
  
  
  completeDomain <- function(x.mat, y.mat, po.x, po.y, min.freq = 1){# ensure sample set covers all time points
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
    
    for (i in 1:n) {
      fp <- fragments()
      po.x[i,fp] <- NA
      po.y[i,fp] <- NA
    }
    
    if(include.full){# include full curve if True
      full.curves <-  sample(1:n, max(1,round((sample(1:20,1)/100)*n)))
      po.x[full.curves,] <- x.mat[full.curves,]
      po.y[full.curves,] <- y.mat[full.curves,]
    }
    
    if(frag.size == "S" & full.domain){
      return(completeDomain(x.mat, y.mat, po.x, po.y))
    }
    
    return(list("po.y"=po.y, "po.x"=po.x))
  }
  
  ################ Base Functions ###################
  
  Kraus <- function(tps){
    K <- 100
    ksi <- sqrt(exp(-(((1:K)-1)^2)/5))*rnorm(K)
    eta <- sqrt(exp(-((1:K)^2)/5))*rnorm(K)
    vk = 3^(-(2*(1:K)-1)); wk = 3^(-(2*(1:K)))
    x <- colSums(t(sapply(1:K, function(k) (sqrt(2)*sqrt(vk[k])*ksi[k]*(cos(2*k*t(pi*tps))/sqrt(5))) + 
                            (sqrt(2)*sqrt(wk[k])*eta[k]*(sin(2*k*t(pi*tps))/sqrt(5))) )))
    return(x)
  }
  
  Alois <- function(tps){
    mu <- tps^(ifelse(base.func[[1]] == 2, 2,1))+sin(2*pi*tps)
    
    K = 50
    ksi.1 <- 50*sqrt(exp(-(((1:K)-1)^2)/ifelse(base.func[[1]] == 2,1,5)))*rnorm(K)
    ksi.2 <- 50*sqrt(exp(-((1:K)^2)/ifelse(base.func[[1]] == 2,1,5)))*rnorm(K)  
    #t((t(ksi.1)%*%(cos((1:K)%*%t(pi*tps))/sqrt(5))) + (t(ksi.2)%*%(sin((1:K)%*%t(pi*tps))/sqrt(5))))
    variations <-colMeans(t(sapply(1:K, function(k) (ksi.1[k]*(cos(k*t(pi*tps))/sqrt(5))) + (ksi.2[k]*(sin(k*t(pi*tps))/sqrt(5))) )))
    
    x <- mu + variations
    if(change.rng) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    return(x)
  }
  
  
  mu.sig <- function(tps){
    kern <- function(s, t) {
      if(is.null(cov.fun)){
        phis <- list(phi.1 = function(t) 1 + t*0, phi.2 = function(t) (2*t - 1)*sqrt(3),
                     phi.3 <- function(t) (6*t^2 - 6*t + 1)*sqrt(5), phi.4 <- function(t) (20 * t^3 - 30*t^2 +12*t -1)*sqrt(7))
        sum(sapply(1:4, function(i) 0.5^(i-1)*phis[[i]](s)*phis[[i]](t)))
      }else if(tolower(cov.fun) =="mercer"){
        MercerKernDecomp(s,t)
      }else if(tolower(cov.fun) == "matern"){
        Matern(s,t)
      }else{
        eval(parse(text = cov.fun))
      }
    }# covariance function
    
    mu <- function (t){
      if(is.null(mean.fun)) t*0
      else eval(parse(text = mean.fun))
    } # mean function
    
    sig <- matrix(nrow = length(tps), ncol = length(tps))
    
    if(is.matrix(cov.fun) & is.numeric(cov.fun)){
      if(ncol(cov.fun) != grid.size) stop("dimension of cov matrix different from grid.size")
      sig <- cov.fun
    }else{
      sig <- sapply(1:length(tps), function(i) 
        sapply(1: length(tps), function(j) kern(tps[i], tps[j])))
    }
    
    if(is.vector(mean.fun) & is.numeric(mean.fun)){
      if(length(mean.fun) != grid.size) stop("dimension of mean vector length different from grid.size")
      m <- mean.fun
    }else{
      m <- mu(tps)
    }
    
    return(list("mean" = m, "cov" = sig))
  }
  
  GP <- function(tps){
    GPparams <- mu.sig(tps)
    x <- MASS::mvrnorm(1, GPparams[["mean"]], GPparams[["cov"]])
    if(change.rng) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    return(x)
  }
  
  
  compositeFunctions <- function(Lt){
    
    if(irreg.domain & POFD.type != "sparse"){
      x <- vector("list", n)
      z <- runif(n, -2/3, 2/3)
      m <- n/5
      x <- lapply(1:n, function(i) {
        if(i <= m) (Lt[[i]]+sin(2*pi*Lt[[i]]))*z[i]
        else if(i > m & i <= m*2) (7*(Lt[[i]]-0.5)^3)*z[i]
        else if(i > m*2 & i <= m*3) (0.5*exp(Lt[[i]])+ (Lt[[i]]-0.5)^2)*z[i]
        else if(i > m*3 & i <= m*3) (Lt[[i]]-cos(2*pi*Lt[[i]]))*z[i]
        else (sin(2*pi*Lt[[i]])+cos(2*pi*Lt[[i]]))*z[i]
      })
      
      if(change.rng) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
      return(x)
    }else{
      t <- Lt
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
      if(change.rng) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
      return(t(x))
    }
    
  }
  
  
  Delaigle <- function(tps, xt = base.func$optn.args$xt, mu = base.func$optn.args$mu, class= 1){# based on simulations in Delaigle and Hall 2013
    
    if(class == 2){
      if(mu == 1) {c <- 12; d <- 4}
      else {c <- 40; d <- 12}
    }else{
      if(mu == 1) {c <- 15; d <- 5}
      else {c <- 50; d <- 20}
    }
    
    if(is.null(base.func$optn.args)) stop("Specify which xt and mu")
    if(mu == 1){
      mu <- function(t) sin(t/c)/(((0.1*t-d)^2)+1)
    }else {
      mu <- function(t) inv.logit((t-c)/d)
    }
    
    if(xt == 1) {
      z <- runif(1, -2/3, 2/3)
      f <- function(t, u) 0.02*((3*t+100)*(u+1))^(1/2)
      x <-  mu(tps) + z*f(tps,z)
    }else if(xt == 2){
      s <- runif(1, -5, 10); u <- runif(1, -1,1)
      v <- runif(1, 0.025, 0.05); w <- runif(1, 2,3)
      z <- runif(1, 0.1, 0.5)
      x <-  mu(tps-s) + (u+v*sin(tps/w))*(z+sin(tps*(10^3)*pi))
    } else{
      u <- runif(1, -1,1); z <- rnorm(1, 0, 0.04)
      v <- runif(1, 0.025, 0.05); w <- runif(1, 2,3)
      x <- mu(tps) +  u+v*sin(t/w)+z
    }
    
    if(change.rng) x <- range.norm(x, a=norm.range[2], b = norm.range[3])
    return (x)
    
  }
  
  
  ########## Covariance Functions ##############
  
  
  MercerKernDecomp <- function(s, t){# finite basis expansion
    if(is.null(args.Mercer)) stop("Specify arguments for Mercer's kernel decomposition")
    if(length(args.Mercer) < 5) stop("Incorrect number of args specified Mercer's kernel decomposition")
    k <- args.Mercer$k; lambda = args.Mercer$lambda
    phi <- args.Mercer$phi; alternate.k = args.Mercer$alternate.k; repeat.phi = args.Mercer$repeat.phi
    if(length(phi) == 1) phi = list(phi)
    
    if(!repeat.phi & length(k) != length(phi)) stop("In args.Mercer, number of basis functions not equal lenght of K and repeat.phi set to FALSE")
    
    if(length(k) == length(phi)){
      phi.list <- vector(mode = "list", length = length(phi))
      for(i in 1:length(k)){
        phi.list[[i]] <- function(t) eval(parse(text = phi[[i]]))
      }
      return(sum(sapply(k, function(i) eval(parse(text = lambda))*phi.list[[i]](s)*phi.list[[i]](t))))
    }
    
    if(length(k) != length(phi) & repeat.phi){
      phi.list <- vector(mode = "list", length = length(phi))
      for(i in 1:length(phi)){
        phi.list[[i]] <- function(t) eval(parse(text = phi[[i]]))
      }
      if(alternate.k){
        if(length(phi) > 2) stop("In args.Mercer, only two phi support for alternating.k at this time")
        phi1 <- sapply(which(k%%2 == 1), function(i) eval(parse(text = lambda))*phi.list[[1]](s)*phi.list[[1]](t))
        phi2 <- sapply(which(k%%2 == 0), function(i) eval(parse(text = lambda))*phi.list[[2]](s)*phi.list[[2]](t))
        return(sum(phi1,phi2))
      }else{
        return(sum(sapply(k, function(i) {
          sum(sapply(1:length(phi), function(j) eval(parse(text = lambda))*phi.list[[j]](s)*phi.list[[j]](t)))
        })))
      }
    }
  }
  
  
  Matern <- function(s,t){
    d <- abs(s-t)
    if(is.null(Matern.args)) {
      l = 1; v = 1
    }else{
      l = Matern.args$l; v = Matern.args$v
    }
    d=ifelse(d == 0,1e-20,d)
    res = (2^(1-v) / gamma(v))*((sqrt(2*v) * d / l)^v)*besselK(sqrt(2*v) * d / l, nu = v)
    return(res)
  }
  
  
  ############## Set with class ############
  
  makeClass <- function(){# generate class
    C1 <- sample(1:n, ceiling(n/2))
    C2 <- setdiff(1:n, C1)
    K <- c(rep(1, length(C1)), rep(2, length(C2)))
    sec1 <- 1:length(C1); sec2 <- (length(C1)+1):n
    
    if(irreg.domain & POFD.type != 'sparse'){
      Lt <- irregFragment()
      Lx1 <- lapply(Lt, function(t) Delaigle(t))[C1]
      Lx2 <- lapply(Lt, function(t) Delaigle(t, class = 2))[C2]
      Lt <- c(Lt[C1], Lt[C2])
      Ly1 <- lapply(Lx1, function(x) x+rnorm(length(x), 0, err.sd))
      Ly2 <- lapply(Lx2, function(x) x+rnorm(length(x), 0, err.sd))
      
      Lx <- c(Lx1, Lx2); Ly <- c(Ly1, Ly2)
      x <- list2mat(Lx)
      y <- list2mat(Ly)
      
      allx <- list("True.Functions" = x, "True.Functions.C1" = x[sec1,], "True.Functions.C2" = x[sec2,])
      #ally <- list("True.Functions.C1" = x[sec1,], "True.Functions.C2" = x[sec2,])
      mu <- list("mu.C1" = colMeans(allx[[2]]), "mu.C2" = colMeans(allx[[3]]))
      
      Ts <- round(seq(grid.range[1],grid.range[2], length.out = grid.size),4)
      Ts.indexes <- lapply(Lt, function(t) which(t %in% Ts))
      Lt.reg <- rep(list(Ts), n)
      Lx.reg <- lapply(1:n, function(i) Lx[[i]][Ts.indexes[[i]]])
      Ly.reg <- lapply(1:n, function(i) Ly[[i]][Ts.indexes[[i]]])
      pofd.x <- list("Lx" = lapply(1:n, function(i) Lx[[i]][-Ts.indexes[[i]]]), "Lt" = lapply(1:n, function(i) Lt[[i]][-Ts.indexes[[i]]]))
      pofd.y <- list("Ly" = lapply(1:n, function(i) Ly[[i]][-Ts.indexes[[i]]]), "Lt" = lapply(1:n, function(i) Lt[[i]][-Ts.indexes[[i]]]))
      all.obs.pts <- sort(unique(unlist(Lt)))
      cov.x = cov(x)
      cov.pd <- as.matrix(Matrix::nearPD(cov.x)$mat)
      if(base.func[[1]] == 4) pop.mu.cov <- mu.sig(Ts)
      else pop.mu.cov = NULL
      return(list("Grid" = list("All.obs.points" = all.obs.pts, "Reg.Grid" = Ts),  "True.Means" = mu, "True.Functions" = allx, "Dense.Functions" = y, 
                  "True.Covariance" = cov.x, "True.Covariance.PD" = cov.pd,"True.List" = list("Lx" = Lx.reg, "Lt" = Lt.reg), "pop.mu.cov" = pop.mu.cov,
                  "Dense.List" = list("Ly" = Ly.reg, "Lt" = Lt.reg),"POFDs" = pofd.y, "POFDs.True.Functions" = pofd.x, "POFD.args" = call.args, "classes" = K))
      
      
    }else{
      Ts <- seq(grid.range[1],grid.range[2], length.out = grid.size)
      err.mat <- matrix(rnorm(n*grid.size,0, err.sd), nrow = n, ncol = grid.size)
      x1 <- sapply(1:n, function(i) Delaigle(Ts))[,C1]
      x2 <- sapply(1:n, function(i) Delaigle(Ts, class = 2))[,C2]
      x1 <- t(x1); x2 <- t(x2)
      y1 <- x1+err.mat[C1,] ; y2 <- x2+err.mat[C2,]
      x <- rbind(x1, x2)
      y <- rbind(y1, y2)
      allx <- list("True.Functions" = x, "True.Functions.C1" = x1, "True.Functions.C2" = x2)
      mu <- list("mu.C1" = colMeans(x1), "mu.C2" = colMeans(x2))
      if(POFD.type == "sparse"){
        pofd <- sparsePOFD(x.mat=x, y.mat=y)
      }else{
        pofd <- fragPOFD(x.mat=x, y.mat=y)
      }
      cov.x = cov(x)
      cov.pd <- as.matrix(Matrix::nearPD(cov.x)$mat)
      if(base.func[[1]] == 4) pop.mu.cov <- mu.sig(Ts)
      else pop.mu.cov = NULL
      return(list("Grid" = Ts, "True.Means" = mu, "True.Functions" = allx, "Dense.Functions" = y, "Dense.List" = mat2List(y),
                  "True.Covariance" = cov.x, "True.Covariance.PD" = cov.pd, "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x,
                  "pop.mu.cov" = pop.mu.cov, "POFD.List" = mat2List(pofd$po.y), "POFD.args" = call.args, "classes" = K))
    }
  }
  
  
  ########### Make sample set ###########
  
  generateSim <- function(){
    
    if(classify) makeClass()
    else{
      if(irreg.domain & POFD.type != 'sparse'){
        Lt <- irregFragment()
        Lx <- switch(base.func[[1]],
                     lapply(Lt, function(t) Alois(t)),
                     lapply(Lt, function(t) Alois(t)),
                     lapply(Lt, function(t) Kraus(t)),
                     lapply(Lt, function(t) GP(t)),
                     compositeFunctions(Lt),
                     lapply(Lt, function(t) Delaigle(t))
        )
        Ly <- lapply(Lx, function(x) x+rnorm(length(x), 0, err.sd))
        x <- list2mat(Lx)
        y <- list2mat(Ly)
        Ts <- round(seq(grid.range[1],grid.range[2], length.out = grid.size),4)
        Ts.indexes <- lapply(Lt, function(t) which(t %in% Ts))
        Lt.reg <- rep(list(Ts), n)
        Lx.reg <- lapply(1:n, function(i) Lx[[i]][Ts.indexes[[i]]])
        Ly.reg <- lapply(1:n, function(i) Ly[[i]][Ts.indexes[[i]]])
        pofd.x <- list("Lx" = lapply(1:n, function(i) Lx[[i]][-Ts.indexes[[i]]]), "Lt" = lapply(1:n, function(i) Lt[[i]][-Ts.indexes[[i]]]))
        pofd.y <- list("Ly" = lapply(1:n, function(i) Ly[[i]][-Ts.indexes[[i]]]), "Lt" = lapply(1:n, function(i) Lt[[i]][-Ts.indexes[[i]]]))
        all.obs.pts <- sort(unique(unlist(Lt)))
        mu <- colMeans(x); cov.x = cov(x)
        cov.pd <- as.matrix(Matrix::nearPD(cov.x)$mat)
        if(base.func[[1]] == 4) pop.mu.cov <- mu.sig(Ts)
        else pop.mu.cov = NULL
        return(list("Grid" = list("All.obs.points" = all.obs.pts, "Reg.Grid" = Ts),  "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, 
                    "True.Covariance" = cov(x), "True.Covariance.PD" = cov.pd, "pop.mu.cov" = pop.mu.cov,
                    "True.List" = list("Lx" = Lx.reg, "Lt" = Lt.reg), "Dense.List" = list("Ly" = Ly.reg, "Lt" = Lt.reg),
                    "POFDs" = pofd.y, "POFDs.True.Functions" = pofd.x, "POFD.args" = call.args))
      }else{
        Ts <- seq(grid.range[1], grid.range[2], length.out = grid.size)
        err.mat <- matrix(rnorm(n*grid.size,0, err.sd), nrow = n, ncol = grid.size)
        x <- switch(base.func[[1]],
                    sapply(1:n, function(i) Alois(Ts)),
                    sapply(1:n, function(i) Alois(Ts)),
                    sapply(1:n, function(i) Kraus(Ts)),
                    sapply(1:n, function(i) GP(Ts)),
                    compositeFunctions(Ts),
                    sapply(1:n, function(i) Delaigle(Ts))
        )
        x <- t(x)
        y <- x+err.mat
        if(POFD.type == "sparse"){
          pofd <- sparsePOFD(x.mat=x, y.mat=y)
        }else{
          pofd <- fragPOFD(x.mat=x, y.mat=y)
        }
        mu <- colMeans(x); cov.x = cov(x)
        cov.pd <- as.matrix(Matrix::nearPD(cov.x)$mat)
        if(base.func[[1]] == 4) pop.mu.cov <- mu.sig(Ts)
        else pop.mu.cov = NULL
        return(list("Grid" = Ts, "True.Mean" = mu, "True.Functions" = x, "Dense.Functions" = y, "Dense.List" = mat2List(y), "True.Covariance" = cov.x, "pop.mu.cov" = pop.mu.cov,
                    "True.Covariance.PD" = cov.pd, "POFDs" = pofd$po.y, "POFDs.True.Functions" = pofd$po.x,"POFD.List" = mat2List(pofd$po.y), "POFD.args" = call.args))
      }
    }
    
  }
  
  sim <- generateSim()
  class(sim) = c("POFD", "list")
  
  return(sim)
  
  
  
}#POFD.sim


