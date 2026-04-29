options(scipen=999)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(brms)
library(DHARMa)
library(DHARMa.helpers)

effectsizes <- data.frame(cs=0.05*((1:100)/100))

r <- readRDS("compnet_results_bin.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    effectsizes[i, "phydist_bin_mean"] <- mean(r[[i]])
    effectsizes[i, "phydist_bin_lower"] <- quantile(r[[i]], 0.025)
    effectsizes[i, "phydist_bin_upper"] <- quantile(r[[i]], 0.975)
  }
  if(class(r[[i]])=="character"){
    effectsizes[i, "phydist_bin_mean"] <- NA
    effectsizes[i, "phydist_bin_lower"] <- NA
    effectsizes[i, "phydist_bin_upper"] <- NA
  }
}

r <- readRDS("compnet_results_fnch.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    effectsizes[i, "phydist_fnch_mean"] <- mean(r[[i]])
    effectsizes[i, "phydist_fnch_lower"] <- quantile(r[[i]], 0.025)
    effectsizes[i, "phydist_fnch_upper"] <- quantile(r[[i]], 0.975)
  }
  if(class(r[[i]])=="character"){
    effectsizes[i, "phydist_fnch_mean"] <- NA
    effectsizes[i, "phydist_fnch_lower"] <- NA
    effectsizes[i, "phydist_fnch_upper"] <- NA
  }
}






set.seed(32190)
# Create alternate 'meta-analysis datasets', in which posterior
# estimates from each stage-1 analysis are thinned for computational
# efficiency in stage 2:
bin_posteriors <- data.frame(cs=c(),
                             draw=c())
bin_posteriors_subsamp <- data.frame(cs=c(),
                                     draw=c())
r <- readRDS("compnet_results_bin.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    tmp <- data.frame(cs=i/2000,
                      draw=r[[i]])
    bin_posteriors <- rbind(bin_posteriors, tmp)
    bin_posteriors_subsamp <- rbind(bin_posteriors_subsamp, tmp[sample(1:nrow(tmp), size=50, replace=F),])
  }
}

fnch_posteriors <- data.frame(cs=c(),
                             draw=c())
fnch_posteriors_subsamp <- data.frame(cs=c(),
                                      draw=c())
r <- readRDS("compnet_results_fnch.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    tmp <- data.frame(cs=i/2000,
                      draw=r[[i]])
    fnch_posteriors <- rbind(fnch_posteriors, tmp)
    fnch_posteriors_subsamp <- rbind(fnch_posteriors_subsamp, tmp[sample(1:nrow(tmp), size=50, replace=F),])
  }
}
# Create factor variables for 'dataset' random effect in meta-analysis
bin_posteriors_subsamp$cs_fac <- factor(bin_posteriors_subsamp$cs)
fnch_posteriors_subsamp$cs_fac <- factor(fnch_posteriors_subsamp$cs)


# GLMM meta-analyses of effect size estimate vs. true competition strength
# Gaussian errors
# Linear effect of true strength on both mean and variance (to account for heteroscedasticity)
# Random effect for 'dataset'
bin_mod <- brm(bf(draw ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                      data=bin_posteriors_subsamp,
                      chains=2,
                      cores=2,
                      warmup=1000,
                      iter=10000,
                      family=student(),
                      seed=2894,
                      control=list(adapt_delta=0.95))
bin_mod_simres <- dh_check_brms(bin_mod, integer=FALSE)
testUniformity(bin_mod_simres)
testQuantiles(bin_mod_simres)
summary(bin_mod)
conditional_effects(bin_mod)
saveRDS(bin_mod, "bin_mod.RDS")
bin_mod <- readRDS("bin_mod.RDS")


fnch_mod <- brm(bf(draw ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                       data=fnch_posteriors_subsamp,
                       chains=2,
                       cores=2,
                       warmup=1000,
                       iter=5000,
                       family=student(),
                       seed=2894,
                       control=list(adapt_delta=0.99))
fnch_mod_simres <- dh_check_brms(fnch_mod, integer=FALSE)
testUniformity(fnch_mod_simres)
testQuantiles(fnch_mod_simres)
summary(fnch_mod)
conditional_effects(fnch_mod)
saveRDS(fnch_mod, "fnch_mod.RDS")
fnch_mod <- readRDS("fnch_mod.RDS")





bin_mod_eff <- conditional_effects(bin_mod)$cs
fnch_mod_eff <- conditional_effects(fnch_mod)$cs

# Save outputs in the other metacommunity simulation study simulation directory
# so the results can all go into one figure
saveRDS(effectsizes, "../compnet_ms_metacom_sim_effectsizes/effectsizes_phy.RDS")
saveRDS(bin_mod_eff, "../compnet_ms_metacom_sim_effectsizes/bin_mod_eff_phy.RDS")
saveRDS(fnch_mod_eff, "../compnet_ms_metacom_sim_effectsizes/fnch_mod_eff_phy.RDS")
saveRDS(bin_posteriors, "../compnet_ms_metacom_sim_effectsizes/bin_posteriors_phy.RDS")
saveRDS(fnch_posteriors, "../compnet_ms_metacom_sim_effectsizes/fnch_posteriors_phy.RDS")
