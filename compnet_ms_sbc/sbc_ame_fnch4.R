### This script simulates data from a rank-2 compnet model
# and then analyzes the data with compnet to confirm that
# compnet can recover correct parameter values.

set.seed(3290)

options(scipen=999) # preserve plenty of numerical precision
library(compnet)
library(rstan)
library(mvtnorm)
library(foreach)
library(doParallel)
registerDoParallel(16)

# Set up data structure
dyads <- t(utils::combn(c(1:20),2))
d <- as.data.frame(dyads)
names(d) <- c("spAid", "spBid")

# For simplicity, assume each species occupies 5 sites, and 100 sites are available.
# Use BiasedUrn notation:
spp <- sort(unique(c(d$spAid, d$spBid)))

# Randomly draw model parameter values
int <- 0.8
sigma <- 1.4
sigma_olre <- 0.8
trait1_trait2_corr <- 0
beta_trait1_sp <- -0.3
beta_trait1_dist <- 1.3
beta_trait2_sp <- -0.1
beta_trait2_dist <- 1.5

# The parameters below are not uniquely identifiable
# due to rotational invariance, reflectional invariance,
# label-switching (when model rank>1), and redundant
# scale factors between lambda_diag and multi_effects_sigma.
# Instead of fitting all these parameters, compnet just z-scores
# each latent factor across species and treats each latent factor
# as independent. My intent is to show that even when we specify all
# these nonidentifiable parameters in the data-generating process, the key quantity they
# produce in the linear predictor--i.e., the product of U, the matrix of latent
# factor values, lambda, and U-transpose--is still identified by my simpler
# model with orthogonal, z-scored latent factors. This approach deals with the 
# rotational invariance issue and the redundant scale factors.
# ame_binomial.stan also fixes the first element of each latent factor to be positive,
# to deal with reflectional invariance, and sorts lambda_diag, to deal with
# label-switching.
multi_effects_corr <- 0
multi_effects_sigma <- c(0.4, 0.4)
multi_effects_sigma <- diag(multi_effects_sigma) %*% 
  matrix(c(1, multi_effects_corr, multi_effects_corr, 1), nrow=2) %*% 
  diag(multi_effects_sigma)
lambda_diag <- c(-1.5, -1.3)


# Simulate one data set for a dummy model (see below)
traits <- rmvnorm(n=length(spp), mean=c(0,0), sigma=matrix(c(1, trait1_trait2_corr, trait1_trait2_corr, 1), nrow=2))
latfacs <- rmvnorm(n=length(spp), mean=c(0,0), sigma=multi_effects_sigma)
spdata <- data.frame(spid=spp,
                     a=rnorm(n=length(spp), mean=0, sd=sigma),
                     latfac1=latfacs[,1],
                     latfac2=latfacs[,2],
                     trait1=scale(traits[,1]),
                     trait2=scale(traits[,2]))
for(i in 1:nrow(d)){
  d[i, "trait1_A"] <- spdata[spdata$spid==d[i,"spAid"], "trait1"]
  d[i, "trait1_B"] <- spdata[spdata$spid==d[i,"spBid"], "trait1"]
  d[i, "trait2_A"] <- spdata[spdata$spid==d[i,"spAid"], "trait2"]
  d[i, "trait2_B"] <- spdata[spdata$spid==d[i,"spBid"], "trait2"]
  d[i, "avalA"] <- spdata[spdata$spid==d[i,"spAid"], "a"]
  d[i, "avalB"] <- spdata[spdata$spid==d[i,"spBid"], "a"]
  d[i, "latfac1A"] <- spdata[spdata$spid==d[i,"spAid"], "latfac1"]
  d[i, "latfac1B"] <- spdata[spdata$spid==d[i,"spBid"], "latfac1"]
  d[i, "latfac2A"] <- spdata[spdata$spid==d[i,"spAid"], "latfac2"]
  d[i, "latfac2B"] <- spdata[spdata$spid==d[i,"spBid"], "latfac2"]
}

for(i in 1:nrow(d)){
  d[i, "latfacterm"] <- (as.matrix(d[i,c("latfac1A", "latfac2A")]))%*%diag(lambda_diag, nrow=2)%*%(t(d[i,c("latfac1B", "latfac2B")]))
}

d$olre <- rnorm(n=nrow(d), sd=sigma_olre)

d$trait1_dist <- abs(d$trait1_A-d$trait1_B)
d$trait2_dist <- abs(d$trait2_A-d$trait2_B)
d$alpha <- int + d$avalA + d$avalB +
  beta_trait1_sp*d$trait1_A +
  beta_trait1_sp*d$trait1_B +
  beta_trait1_dist*d$trait1_dist+
  beta_trait2_sp*d$trait2_A +
  beta_trait2_sp*d$trait2_B +
  beta_trait2_dist*d$trait2_dist+
  d$latfacterm+
  d$olre

d$both <- sapply(exp(d$alpha), BiasedUrn::rFNCHypergeo, nran=1, m1=5, m2=95, n=5)

# Set up dummy model so models can be fit without recompiling in the foreach loop
datalist <- list(
  n_nodes=length(spp),
  sp_occ=rep(5, times=length(spp)),
  n_sites=100,
  N=nrow(d),
  Xdy_cols=2,
  Xsp_cols=2,
  spAid=d$spAid,
  spBid=d$spBid,
  Xdy=d[c("trait1_dist", "trait2_dist")],
  XA=d[c("trait1_A", "trait2_A")],
  XB=d[c("trait1_B", "trait2_B")],
  K=2,
  both=d$both,
  either=d$either,
  prior_intercept_scale=5,
  prior_betas_scale=5,
  prior_sigma_addeff_rate=1,
  prior_lambda_scale=5,
  prior_sigma_olre_rate=1)
stanmod_dummy <- stan(file="ame_fnchypg.stan",
                      data=datalist,
                      cores=1,
                      chains=1,
                      warmup=1,
                      iter=2,
                      verbose=F)



# Simulate and analyze 16 replicate data sets
sbc_mods_and_data <- foreach(u=1:16) %dopar% {
  library(compnet)
  library(rstan)
  library(mvtnorm)
  dyads <- t(utils::combn(c(1:20),2))
  d <- as.data.frame(dyads)
  names(d) <- c("spAid", "spBid")
  
  d$either <- 100
  spp <- sort(unique(c(d$spAid, d$spBid)))
  
  traits <- rmvnorm(n=length(spp), mean=c(0,0), sigma=matrix(c(1, trait1_trait2_corr, trait1_trait2_corr, 1), nrow=2))
  latfacs <- rmvnorm(n=length(spp), mean=c(0,0), sigma=multi_effects_sigma)
  spdata <- data.frame(spid=spp,
                       a=rnorm(n=length(spp), mean=0, sd=sigma),
                       latfac1=latfacs[,1],
                       latfac2=latfacs[,2],
                       trait1=scale(traits[,1]),
                       trait2=scale(traits[,2]))
  for(i in 1:nrow(d)){
    d[i, "trait1_A"] <- spdata[spdata$spid==d[i,"spAid"], "trait1"]
    d[i, "trait1_B"] <- spdata[spdata$spid==d[i,"spBid"], "trait1"]
    d[i, "trait2_A"] <- spdata[spdata$spid==d[i,"spAid"], "trait2"]
    d[i, "trait2_B"] <- spdata[spdata$spid==d[i,"spBid"], "trait2"]
    d[i, "avalA"] <- spdata[spdata$spid==d[i,"spAid"], "a"]
    d[i, "avalB"] <- spdata[spdata$spid==d[i,"spBid"], "a"]
    d[i, "latfac1A"] <- spdata[spdata$spid==d[i,"spAid"], "latfac1"]
    d[i, "latfac1B"] <- spdata[spdata$spid==d[i,"spBid"], "latfac1"]
    d[i, "latfac2A"] <- spdata[spdata$spid==d[i,"spAid"], "latfac2"]
    d[i, "latfac2B"] <- spdata[spdata$spid==d[i,"spBid"], "latfac2"]
  }
  
  for(i in 1:nrow(d)){
    d[i, "latfacterm"] <- (as.matrix(d[i,c("latfac1A", "latfac2A")]))%*%diag(lambda_diag, nrow=2)%*%(t(d[i,c("latfac1B", "latfac2B")]))
  }
  
  d$olre <- rnorm(n=nrow(d), sd=sigma_olre)
  
  d$trait1_dist <- abs(d$trait1_A-d$trait1_B)
  d$trait2_dist <- abs(d$trait2_A-d$trait2_B)
  d$alpha <- int + d$avalA + d$avalB +
    beta_trait1_sp*d$trait1_A +
    beta_trait1_sp*d$trait1_B +
    beta_trait1_dist*d$trait1_dist+
    beta_trait2_sp*d$trait2_A +
    beta_trait2_sp*d$trait2_B +
    beta_trait2_dist*d$trait2_dist+
    d$latfacterm+
    d$olre
  
  d$both <- sapply(exp(d$alpha), BiasedUrn::rFNCHypergeo, nran=1, m1=5, m2=95, n=5)
  
  datalist <- list(
    n_nodes=length(spp),
    sp_occ=rep(5, times=length(spp)),
    n_sites=100,
    N=nrow(d),
    Xdy_cols=2,
    Xsp_cols=2,
    spAid=d$spAid,
    spBid=d$spBid,
    Xdy=d[c("trait1_dist", "trait2_dist")],
    XA=d[c("trait1_A", "trait2_A")],
    XB=d[c("trait1_B", "trait2_B")],
    K=2,
    both=d$both,
    either=d$either,
    prior_intercept_scale=5,
    prior_betas_scale=5,
    prior_sigma_addeff_rate=1,
    prior_lambda_scale=5,
    prior_sigma_olre_rate=1)
  
  stanmod <- stan(file="ame_fnchypg.stan",
                  fit=stanmod_dummy,
                  pars=c("intercept", # Marginalize over the individual values of each latent factor ("U"), since these are known to be nonidentifiable due to reflectional and rotational invariance, as well as potential label-switching among multiple latent factor dimensions when K>1. The Stan code has tools designed to deal with this, but only for computational efficiency reasons; it's not a problem for inference.
                         "beta_dy",
                         "beta_sp",
                         "sigma",
                         "a",
                         "sigma_olre",
                         "olre",
                         "alpha",
                         "latfacterm", # In contrast to individual species' U values, ULU values (sensu Hoff and colleagues) need to be identifiable, so we'll include them to make sure they sample well.
                         "lp__"),
                  data=datalist,
                  cores=1,
                  chains=1,
                  warmup=1000,
                  iter=3000,
                  verbose=F,
                  control=list(adapt_delta=0.99))
  truepars <- list(int, sigma, sigma_olre, beta_trait1_sp, beta_trait1_dist, beta_trait2_sp, beta_trait2_dist, latfacs)
  tmpout <- list(stanmod, d, datalist, truepars)
  tmpout
}
saveRDS(sbc_mods_and_data, "sbc_ame_fnch4_out.RDS")
