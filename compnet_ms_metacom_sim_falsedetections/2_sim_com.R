options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(10) # adjust number of cores for different machine

# Set the number of species
species <- 100
saveRDS(species, "species.RDS")

# Set the maximum population growth rate, rmax
rmax <- 5

set.seed(490)
res <- foreach(u=1:100) %dopar% {
  i <- sample(1:30, size=1)
  
  # function for density-independent growth rate of species i at site x at time t
  # rmax is a single value.
  # sp_df is a data frame, in which each species is a row.
  # sp_df contains vectors for the "trait value" zi (i.e., optimal environment),
  # as well as niche breadth sigmai.
  # env_df_t is a data frame, in which each site is a row.
  # env_df_t contains a vector for the environmental condition at each site at time t.
  # produces a matrix, in which rows are sites and columns are species.
  COMP_r_t0 <- function(rmax, sp_df, env_df_t){
    df <- merge(sp_df, env_df_t)
    df$rix_t0 <- rmax*exp(-((df$zi-df$envx_t)/(2*df$sigmai))^2)
    df <- df[order(rownames(df)),]
    df <- df[order(df$sp),]
    matrix(df$rix_t0, nrow=nrow(env_df_t), ncol=nrow(sp_df))
  }
  
  # function for abundance of species i at site x at time t+1
  # before accounting for immigration or emigration.
  # N_t0 is a matrix in which rows are sites and columns are species.
  # r_t0_mat is a matrix in which rows are sites and columns are species.
  # alpha_mat is a matrix containing the competitive effect of the
  # row species on the column species.
  # produces a matrix in which rows are sites and columns are species.
  COMP_Nhat_t1 <- function(N_t0, r_t0_mat, alpha_mat){
    denom <- N_t0%*%alpha_mat+1
    lambda <- N_t0*(r_t0_mat/denom)
    lambda[lambda<0] <- 0
    vals <- rpois(n=length(lambda), lambda=lambda)
    matrix(vals, nrow=nrow(N_t0), ncol=ncol(N_t0))
  }
  
  # function for number of emigrating individuals of species i from site x at time t
  # Nhat_t1 is a matrix of pre-dispersal abundances at time t1,
  # in which rows are sites and columns are species.
  # a_vec is a vector containing an a value (which describes dispersal propensity) for each species.
  # produces a matrix, in which rows are sites and columns are species.
  COMP_E_t <- function(Nhat_t1, a_vec){
    a_mat <- t(matrix(a_vec))
    a_mat <- a_mat[rep(1, times=nrow(Nhat_t1)),]
    counts <- rbinom(n=length(Nhat_t1), size=Nhat_t1, prob=a_mat)
    matrix(counts, nrow=nrow(Nhat_t1), ncol=ncol(Nhat_t1))
  }
  
  # function for number of immigrating individuals of species i to site x at time t.
  # disp_mat is a dispersal matrix among all pairs of sites (with zero diag).
  # E_t_mat is a matrix of emigrant counts at time t,
  # in which rows are sites and columns are species.
  # produces a matrix, in which rows are sites and columns are species.
  COMP_I_t <- function(disp_mat, E_t){
    round(disp_mat %*% E_t)
  }
  
  # function for full abundance of species i at site x at time t+1
  # all inputs are matrices. rows are sites, columns are species.
  # produces a matrix, in which rows are sites and columns are species.
  COMP_N_t1 <- function(Nhat_t1, E_t, I_t){
    N_t1 <- Nhat_t1 - E_t + I_t
    N_t1[N_t1<0] <- 0
    N_t1
  }
  
  env <- readRDS(paste("thompson_env_", i, ".RDS", sep=""))
  disp_mat <- readRDS(paste("thompson_disp_mat_", i, ".RDS", sep=""))
  
  # Generate data frame of species and trait values
  library(geiger)
  library(phytools)
  sp_df <- data.frame(sp=paste("sp", 1:species, sep="_"))
  tree <- drop.extinct(sim.bdtree(b=2, d=1, stop=c("taxa"), n=species, extinct=FALSE, seed=i))
  tree$tip.label <- sp_df$sp
  # set correlation among traits:
  traitcor <- 0.8
  vcvmat <- matrix(c(1,traitcor,0,0,
                     traitcor,1,0,0,
                     0,0,1,0,
                     0,0,0,1), nrow=4)
  traits <- as.data.frame(sim.char(tree, par=vcvmat, model="BM"))
  names(traits) <- c("dummytrait", "domtrait", "zi", "ndtrait")
  
  # z-scale each trait
  for(j in 1:ncol(traits)){
    traits[,j] <- as.vector(scale(traits[,j], center=TRUE))
  }
  sp_df <- cbind(sp_df, traits)
  
  sp_df$ai <- 0.01
  sp_df$sigmai <- 1
  
  # Generate a matrix of alpha values defining the strength of competition within and between species
  ndtrait_dist <- as.matrix(dist(sp_df$ndtrait))
  domtrait_dist <- t(outer(sp_df$domtrait, sp_df$domtrait, `-`))
  comp_dist <- ndtrait_dist+domtrait_dist
  
  alpha_mat <- 0.05*exp(-comp_dist)
  
  ##########
  N <- list()
  N_init <- matrix(rpois(lambda=0.5, n=nrow(sp_df)*nrow(unique(env[c("x","y")]))), nrow=nrow(unique(env[c("x","y")])), ncol=nrow(sp_df))
  N[[1]] <- N_init
  for(j in 2:200){
    N_t0 <- N[[j-1]]
    if(j/10==round(j/10)){
      N_t0[N_t0==0] <- rpois(lambda=0.5, n=length(N_t0[N_t0==0]))
    }
    env_df_t <- subset(env, time==j)
    r_t0_mat <- COMP_r_t0(rmax=rmax, sp_df=sp_df, env_df_t=env_df_t)
    Nhat_t1 <- COMP_Nhat_t1(N_t0=N_t0, r_t0_mat=r_t0_mat, alpha_mat=alpha_mat)
    E_t <- COMP_E_t(Nhat_t1=Nhat_t1, a_vec=sp_df$ai)
    I_t <- COMP_I_t(disp_mat=disp_mat, E_t=E_t)
    N_t1 <- COMP_N_t1(Nhat_t1=Nhat_t1, E_t=E_t, I_t=I_t)
    N[[j]] <- N_t1
  }
  
  for(j in 201:2200){
    N_t0 <- N[[j-1]]
    env_df_t <- subset(env, time==j)
    r_t0_mat <- COMP_r_t0(rmax=rmax, sp_df=sp_df, env_df_t=env_df_t)
    Nhat_t1 <- COMP_Nhat_t1(N_t0=N_t0, r_t0_mat=r_t0_mat, alpha_mat=alpha_mat)
    E_t <- COMP_E_t(Nhat_t1=Nhat_t1, a_vec=sp_df$ai)
    I_t <- COMP_I_t(disp_mat=disp_mat, E_t=E_t)
    N_t1 <- COMP_N_t1(Nhat_t1=Nhat_t1, E_t=E_t, I_t=I_t)
    N[[j]] <- N_t1
  }
  list(N[[2200]], sp_df)
}
saveRDS(res, "sim_com_output.RDS")
