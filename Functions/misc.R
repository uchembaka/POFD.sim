
POFD.sim2List <- function(FD.mat){
  grid <- seq(0,1, len=ncol(FD.mat)) 
  Ly <- split(t(FD.mat), rep(1:ncol(t(FD.mat)), each = nrow(t(FD.mat))))
  Lt <- lapply(Ly, function(i) grid[!is.na(i)])
  Ly <- lapply(Ly, function(i) i[!is.na(i)])
  return(list(Ly, Lt))
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
    if(length(phi) > 2) stop("Number of basis function must be 2, e.g sin(kt) and cos(kt)")
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