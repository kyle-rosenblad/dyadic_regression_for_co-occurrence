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
    effectsizes[i, "ndtrait_bin_mean"] <- mean(r[[i]][[1]])
    effectsizes[i, "ndtrait_bin_lower"] <- quantile(r[[i]][[1]], 0.025)
    effectsizes[i, "ndtrait_bin_upper"] <- quantile(r[[i]][[1]], 0.975)
    effectsizes[i, "domtrait_bin_mean"] <- mean(r[[i]][[2]])
    effectsizes[i, "domtrait_bin_lower"] <- quantile(r[[i]][[2]], 0.025)
    effectsizes[i, "domtrait_bin_upper"] <- quantile(r[[i]][[2]], 0.975)
  }
  if(class(r[[i]])=="character"){
    effectsizes[i, "ndtrait_bin_mean"] <- NA
    effectsizes[i, "ndtrait_bin_lower"] <- NA
    effectsizes[i, "ndtrait_bin_upper"] <- NA
    effectsizes[i, "domtrait_bin_mean"] <- NA
    effectsizes[i, "domtrait_bin_lower"] <- NA
    effectsizes[i, "domtrait_bin_upper"] <- NA
  }
}


r <- readRDS("compnet_results_fnch.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    effectsizes[i, "ndtrait_fnch_mean"] <- mean(r[[i]][[1]])
    effectsizes[i, "ndtrait_fnch_lower"] <- quantile(r[[i]][[1]], 0.025)
    effectsizes[i, "ndtrait_fnch_upper"] <- quantile(r[[i]][[1]], 0.975)
    effectsizes[i, "domtrait_fnch_mean"] <- mean(r[[i]][[2]])
    effectsizes[i, "domtrait_fnch_lower"] <- quantile(r[[i]][[2]], 0.025)
    effectsizes[i, "domtrait_fnch_upper"] <- quantile(r[[i]][[2]], 0.975)
  }
  if(class(r[[i]])=="character"){
    effectsizes[i, "ndtrait_fnch_mean"] <- NA
    effectsizes[i, "ndtrait_fnch_lower"] <- NA
    effectsizes[i, "ndtrait_fnch_upper"] <- NA
    effectsizes[i, "domtrait_fnch_mean"] <- NA
    effectsizes[i, "domtrait_fnch_lower"] <- NA
    effectsizes[i, "domtrait_fnch_upper"] <- NA
  }
}






set.seed(32190)
# Create alternate 'meta-analysis datasets', in which posterior
# estimates from each stage-1 analysis are thinned for computational
# efficiency in stage 2:
bin_posteriors <- data.frame(cs=c(),
                             draw1=c(),
                             draw2=c())
bin_posteriors_subsamp <- data.frame(cs=c(),
                                     draw1=c(),
                                     draw2=c())
r <- readRDS("compnet_results_bin.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    tmp <- data.frame(cs=i/2000,
                      draw1=r[[i]][[1]],
                      draw2=r[[i]][[2]])
    bin_posteriors <- rbind(bin_posteriors, tmp)
    bin_posteriors_subsamp <- rbind(bin_posteriors_subsamp, tmp[sample(1:nrow(tmp), size=50, replace=F),])
  }
}
fnch_posteriors <- data.frame(cs=c(),
                             draw1=c(),
                             draw2=c())
fnch_posteriors_subsamp <- data.frame(cs=c(),
                                      draw1=c(),
                                      draw2=c())
r <- readRDS("compnet_results_fnch.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    tmp <- data.frame(cs=i/2000,
                      draw1=r[[i]][[1]],
                      draw2=r[[i]][[2]])
    fnch_posteriors <- rbind(fnch_posteriors, tmp)
    fnch_posteriors_subsamp <- rbind(fnch_posteriors_subsamp, tmp[sample(1:nrow(tmp), size=50, replace=F),])
  }
}
# Create factor variables for 'dataset' random effect in meta-analysis
bin_posteriors_subsamp$cs_fac <- factor(bin_posteriors_subsamp$cs)
fnch_posteriors_subsamp$cs_fac <- factor(fnch_posteriors_subsamp$cs)


# GLMM meta-analyses of effect size estimate vs. true competition strength
# Student's errors to account for heavy tails
# Linear effect of true strength on both mean and variance (to account for heteroscedasticity)
# Random effect for 'dataset'
bin_mod1 <- brm(bf(draw1 ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                      data=bin_posteriors_subsamp,
                      chains=2,
                      cores=2,
                      warmup=1000,
                      iter=5000,
                      family=student(),
                      seed=2894,
                      control=list(adapt_delta=0.99))
bin_mod1_simres <- dh_check_brms(bin_mod1, integer=FALSE)
testUniformity(bin_mod1_simres)
testQuantiles(bin_mod1_simres)
summary(bin_mod1)
conditional_effects(bin_mod1)
saveRDS(bin_mod1, "bin_mod1.RDS")
bin_mod1 <- readRDS("bin_mod1.RDS") #


bin_mod2 <- brm(bf(draw2 ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                data=bin_posteriors_subsamp,
                chains=2,
                cores=2,
                warmup=1000,
                iter=4000,
                family=student(),
                seed=2894,
                control=list(adapt_delta=0.99))
bin_mod2_simres <- dh_check_brms(bin_mod2, integer=FALSE)
testUniformity(bin_mod2_simres)
testQuantiles(bin_mod2_simres)
summary(bin_mod2)
conditional_effects(bin_mod2)
saveRDS(bin_mod2, "bin_mod2.RDS")
bin_mod2 <- readRDS("bin_mod2.RDS") #






fnch_mod1 <- brm(bf(draw1 ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                       data=fnch_posteriors_subsamp,
                       chains=2,
                       cores=2,
                       warmup=1000,
                       iter=5000,
                       family=student(),
                       seed=2894,
                       control=list(adapt_delta=0.99))
fnch_mod1_simres <- dh_check_brms(fnch_mod1, integer=FALSE)
testUniformity(fnch_mod1_simres)
testQuantiles(fnch_mod1_simres)
summary(fnch_mod1)
conditional_effects(fnch_mod1)
saveRDS(fnch_mod1, "fnch_mod1.RDS")
fnch_mod1 <- readRDS("fnch_mod1.RDS") #


fnch_mod2 <- brm(bf(draw2 ~ scale(cs) + (1|cs_fac), sigma~s(scale(cs))),
                 data=fnch_posteriors_subsamp,
                 chains=2,
                 cores=2,
                 warmup=1000,
                 iter=5000,
                 family=student(),
                 seed=2894,
                 control=list(adapt_delta=0.99))
fnch_mod2_simres <- dh_check_brms(fnch_mod2, integer=FALSE)
testUniformity(fnch_mod2_simres)
testQuantiles(fnch_mod2_simres)
summary(fnch_mod2)
conditional_effects(fnch_mod2)
saveRDS(fnch_mod2, "fnch_mod2.RDS")
fnch_mod2 <- readRDS("fnch_mod2.RDS") #



bin_mod1_eff <- conditional_effects(bin_mod1)$cs
fnch_mod1_eff <- conditional_effects(fnch_mod1)$cs
bin_mod2_eff <- conditional_effects(bin_mod2)$cs
fnch_mod2_eff <- conditional_effects(fnch_mod2)$cs



# Estimate distribution of Bayesian R2 for each model
bin_mod1_R2 <- bayes_R2(bin_mod1)
fnch_mod1_R2 <- bayes_R2(fnch_mod1)
bin_mod1_R2
fnch_mod1_R2

bin_mod2_R2 <- bayes_R2(bin_mod2)
fnch_mod2_R2 <- bayes_R2(fnch_mod2)
bin_mod2_R2
fnch_mod2_R2



# Estimate posterior distribution of Spearman correlation between competition strength
# and effect size estimate for each model
n_post <- nrow(bin_posteriors) / nrow(effectsizes)
spearman_draws <- data.frame(bin1=c(),
                             fnch1=c(),
                             bin2=c(),
                             fnch2=c())
for(i in 1:n_post){
  cs_draws <- bin_posteriors[i + n_post*(c(0:(nrow(effectsizes)-1))), "cs"]
  bin_effect1_draws <- bin_posteriors[i + n_post*(c(0:(nrow(effectsizes)-1))), "draw1"]
  bin_effect2_draws <- bin_posteriors[i + n_post*(c(0:(nrow(effectsizes)-1))), "draw2"]
  fnch_effect1_draws <- fnch_posteriors[i + n_post*(c(0:(nrow(effectsizes)-1))), "draw1"]
  fnch_effect2_draws <- fnch_posteriors[i + n_post*(c(0:(nrow(effectsizes)-1))), "draw2"]
  
  spearman_draws[i, "bin1"] <- cor(cs_draws, bin_effect1_draws)
  spearman_draws[i, "fnch1"] <- cor(cs_draws, fnch_effect1_draws)
  spearman_draws[i, "bin2"] <- cor(cs_draws, bin_effect2_draws)
  spearman_draws[i, "fnch2"] <- cor(cs_draws, fnch_effect2_draws)
}
head(spearman_draws)

spearman_summ_bin1 <- quantile(spearman_draws$bin1, c(0.025, 0.975))
spearman_summ_fnch1 <- quantile(spearman_draws$fnch1, c(0.025, 0.975))
spearman_summ_bin2 <- quantile(spearman_draws$bin2, c(0.025, 0.975))
spearman_summ_fnch2 <- quantile(spearman_draws$fnch2, c(0.025, 0.975))

spearman_summ_bin1
spearman_summ_fnch1

spearman_summ_bin2
spearman_summ_fnch2

p1 <- ggplot()+
  ggtitle("Binomial Models")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=bin_mod1_eff, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=bin_mod1_eff, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes, aes(x=cs, xend=cs, y=ndtrait_bin_lower, yend=ndtrait_bin_upper),
               linewidth=0.1,
               alpha=0.5)+
  annotate("text", x = 0.001, y = -1.2, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(bin_mod1_R2[3], 2),
                         ", ",
                         round(bin_mod1_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_bin1[1], 2),
                         ", ",
                         round(spearman_summ_bin1[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Niche Differentiation Trait Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title.x=element_blank())
p1

p2 <- ggplot()+
  ggtitle("FNCH Models")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=fnch_mod1_eff, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=fnch_mod1_eff, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes, aes(x=cs, xend=cs, y=ndtrait_fnch_lower, yend=ndtrait_fnch_upper),
               linewidth=0.1, alpha=0.5)+
  annotate("text", x = 0.001, y = -2.6, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(fnch_mod1_R2[3], 2),
                         ", ",
                         round(fnch_mod1_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_fnch1[1], 2),
                         ", ",
                         round(spearman_summ_fnch1[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Niche Differentiation Trait Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title=element_blank())
p2

p3 <- ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=bin_mod2_eff, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=bin_mod2_eff, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes, aes(x=cs, xend=cs, y=domtrait_bin_lower, yend=domtrait_bin_upper),
               linewidth=0.1, alpha=0.5)+
  annotate("text", x = 0.0013, y = -2.1, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(bin_mod2_R2[3], 2),
                         ", ",
                         round(bin_mod2_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_bin2[1], 2),
                         ", ",
                         round(spearman_summ_bin2[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Dominance Trait Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title.x=element_blank())
p3

p4 <- ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=fnch_mod2_eff, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=fnch_mod2_eff, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes, aes(x=cs, xend=cs, y=domtrait_fnch_lower, yend=domtrait_fnch_upper),
               linewidth=0.1, alpha=0.5)+
  annotate("text", x = 0.001, y = -2.45, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(fnch_mod2_R2[3], 2),
                         ", ",
                         round(fnch_mod2_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_fnch2[1], 2),
                         ", ",
                         round(spearman_summ_fnch2[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Dominance Trait Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title=element_blank())
p4


effectsizes_phy <- readRDS("effectsizes_phy.RDS")
bin_mod_eff_phy <- readRDS( "bin_mod_eff_phy.RDS")
fnch_mod_eff_phy <- readRDS( "fnch_mod_eff_phy.RDS")
bin_posteriors_phy <- readRDS( "bin_posteriors_phy.RDS")
fnch_posteriors_phy <- readRDS( "fnch_posteriors_phy.RDS")
bin_mod_phy <- readRDS("../compnet_ms_metacom_sim_effectsizes_phy/bin_mod.RDS")
fnch_mod_phy <- readRDS("../compnet_ms_metacom_sim_effectsizes_phy/fnch_mod.RDS")

bin_mod_phy_R2 <- bayes_R2(bin_mod_phy)
fnch_mod_phy_R2 <- bayes_R2(fnch_mod_phy)

n_post <- nrow(bin_posteriors_phy) / nrow(effectsizes_phy)
spearman_draws_phy <- data.frame(bin=c(),
                             fnch=c())
for(i in 1:n_post){
  cs_draws <- bin_posteriors_phy[i + n_post*(c(0:(nrow(effectsizes_phy)-1))), "cs"]
  bin_effect_phy_draws <- bin_posteriors_phy[i + n_post*(c(0:(nrow(effectsizes_phy)-1))), "draw"]
  fnch_effect_phy_draws <- fnch_posteriors_phy[i + n_post*(c(0:(nrow(effectsizes_phy)-1))), "draw"]

  spearman_draws_phy[i, "bin"] <- cor(cs_draws, bin_effect_phy_draws)
  spearman_draws_phy[i, "fnch"] <- cor(cs_draws, fnch_effect_phy_draws)
}
head(spearman_draws_phy)

spearman_summ_bin_phy <- quantile(spearman_draws_phy$bin, c(0.025, 0.975))
spearman_summ_fnch_phy <- quantile(spearman_draws_phy$fnch, c(0.025, 0.975))


p5 <- ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=bin_mod_eff_phy, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=bin_mod_eff_phy, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes_phy, aes(x=cs, xend=cs, y=phydist_bin_lower, yend=phydist_bin_upper),
               linewidth=0.1, alpha=0.5)+
  annotate("text", x = 0.001, y = 0.365, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(bin_mod_phy_R2[3], 2),
                         ", ",
                         round(bin_mod_phy_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_bin_phy[1], 2),
                         ", ",
                         round(spearman_summ_bin_phy[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Phylogenetic Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1)
p5

p6 <- ggplot()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=fnch_mod_eff_phy, aes(x=cs, ymin=lower__, ymax=upper__), fill="red", alpha=0.25)+
  geom_line(data=fnch_mod_eff_phy, aes(x=cs, y=estimate__), color="red")+
  geom_segment(data=effectsizes_phy, aes(x=cs, xend=cs, y=phydist_fnch_lower, yend=phydist_fnch_upper),
               linewidth=0.1, alpha=0.5)+
  annotate("text", x = 0.002, y = 2.25, size=2.5, hjust=0,
           label = paste("Bayes R²: [",
                         round(fnch_mod_phy_R2[3], 2),
                         ", ",
                         round(fnch_mod_phy_R2[4], 2),
                         "]",
                         "\nRank Corr: [",
                         round(spearman_summ_fnch_phy[1], 2),
                         ", ",
                         round(spearman_summ_fnch_phy[2], 2),
                         "]",
                         sep=""))+
  xlab("Strength of Simulated Competition")+
  ylab("Phylogenetic Distance Effect")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title.y=element_blank())
p6





p1+p2+p3+p4+p5+p6+
  plot_annotation(title="Standardized Effect Estimates")+
  plot_layout(nrow=3)
ggsave("../output_files/fig3.pdf", height=9, width=6)
ggsave("../output_files/fig3.png", height=9, width=6)

