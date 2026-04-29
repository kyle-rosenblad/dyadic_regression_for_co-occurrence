# This script runs the metacommunity simulations described in Thompson et al. (2020) Ecol. Lett.
# I am rebuilding from scratch using the equations described in the paper.
options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(3) # adjust number of cores for different machine

set.seed(12078)
foreach(j=1:30) %dopar% {
  library(sdmTMB)
  library(INLA)
  # Generate spatiotemporal landscape with autocorrelated environmental variable:
  landscape <- unique(round(data.frame(x=runif(100, min = 1, max = 100), y=runif(100, min = 1, max = 100))))
  while(nrow(landscape)<100){
    landscape <- unique(round(data.frame(x=runif(100, min = 1, max = 100), y=runif(100, min = 1, max = 100))))
  }
  time <- data.frame(time=0:2000)
  landscape <- merge(landscape, time)
  mesh <- make_mesh(landscape, xy_cols = c("x", "y"), cutoff = 1)
  env <- sdmTMB_simulate(
    formula = ~ 1,
    data = landscape,
    time="time",
    mesh = mesh,
    family = gaussian(link = "identity"),
    phi=1e-6,
    rho=0.9,
    range = 30,
    sigma_O = 0.5,
    sigma_E = 0.02,
    seed = 34,
    B = 0
  )
  env <- env[c("time", "x", "y", "eta")]
  env$envx_t <- env$eta
  env$eta <- NULL
  rownames(env) <- paste(env$x, env$y, env$time, sep="_")
  
  env$time <- env$time+200
  for(i in 1:199){
    tmp <- subset(env, time==200)
    tmp$time <- i
    rownames(tmp) <- paste(tmp$x, tmp$y, tmp$time, sep="_")
    env <- rbind(tmp, env)
  }
  
  # scale and center environmental variable:
  env$envx_t <- scale(env$envx_t, center=TRUE)
  
  # Generate dispersal matrix:
  dist_mat <- as.matrix(dist(env[env$time==1, c("x", "y")]))
  disp_mat <- exp(-0.1*dist_mat)
  diag(disp_mat) <- 0
  disp_mat <- apply(disp_mat, 1, function(x) x/sum(x))
  
  saveRDS(env, paste("thompson_env_", j, ".RDS", sep=""))
  saveRDS(disp_mat, paste("thompson_disp_mat_", j, ".RDS", sep=""))
}
rm(list=ls())
gc()
