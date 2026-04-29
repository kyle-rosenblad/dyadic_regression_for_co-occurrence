options(scipen=999) # preserve plenty of numerical precision
library(compnet)
library(rstan)
library(foreach)
library(doParallel)
registerDoParallel(cores=8)


library(rstan)
md <- readRDS("sbc_ame_fnch1_out.RDS")
intercept <- md[[1]][[4]][[1]]
sigma <- md[[1]][[4]][[2]]
sigma_olre <- md[[1]][[4]][[3]]
beta_trait1_sp <- md[[1]][[4]][[4]]
beta_trait1_dist <- md[[1]][[4]][[5]]
beta_trait2_sp <- md[[1]][[4]][[6]]
beta_trait2_dist <- md[[1]][[4]][[7]]
for(v in 1:16){
  srmpars <- data.frame(rep=c(),
                        intercept=c(),
                        sigma=c(),
                        sigma_olre=c(),
                        beta_trait1_sp=c(),
                        beta_trait1_dist=c(),
                        beta_trait2_sp=c(),
                        beta_trait2_dist=c())
  for(i in 1:length(md)){
    tmpmod <- md[[i]][[1]]
    tmpdraws <- extract(tmpmod)
    srmpars_tmp <- data.frame(rep=i,
                              intercept=tmpdraws$intercept,
                              sigma=tmpdraws$sigma,
                              sigma_olre=tmpdraws$sigma_olre,
                              beta_trait1_sp=tmpdraws$beta_sp[,1],
                              beta_trait1_dist=tmpdraws$beta_dy[,1],
                              beta_trait2_sp=tmpdraws$beta_sp[,2],
                              beta_trait2_dist=tmpdraws$beta_dy[,2])
    srmpars <- rbind(srmpars, srmpars_tmp)
  }
}
library(ggplot2)
intercept_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=intercept, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=intercept), color="blue", linewidth=1.5)+
  geom_vline(xintercept=intercept)+
  geom_vline(xintercept=mean(srmpars$intercept), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Intercept")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

sigma_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=sigma, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=sigma), color="blue", linewidth=1.5)+
  geom_vline(xintercept=sigma)+
  geom_vline(xintercept=mean(srmpars$sigma), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Additive Species Effect SD")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

sigma_olre_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=sigma_olre, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=sigma_olre), color="blue", linewidth=1.5)+
  geom_vline(xintercept=sigma_olre)+
  geom_vline(xintercept=mean(srmpars$sigma_olre), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Observation-Level Random Effect SD")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

beta_trait1_sp_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=beta_trait1_sp, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=beta_trait1_sp), color="blue", linewidth=1.5)+
  geom_vline(xintercept=beta_trait1_sp)+
  geom_vline(xintercept=mean(srmpars$beta_trait1_sp), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Species Effect, Trait 1")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

beta_trait1_dist_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=beta_trait1_dist, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=beta_trait1_dist), color="blue", linewidth=1.5)+
  geom_vline(xintercept=beta_trait1_dist)+
  geom_vline(xintercept=mean(srmpars$beta_trait1_dist), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Distance Effect, Trait 1")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

beta_trait2_sp_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=beta_trait2_sp, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=beta_trait2_sp), color="blue", linewidth=1.5)+
  geom_vline(xintercept=beta_trait2_sp)+
  geom_vline(xintercept=mean(srmpars$beta_trait2_sp), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Species Effect, Trait 2")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

beta_trait2_dist_plot <- ggplot(srmpars)+
  geom_hline(yintercept=0)+
  geom_density(aes(x=beta_trait2_dist, group=rep), color=alpha("blue", 0.1))+
  geom_density(aes(x=beta_trait2_dist), color="blue", linewidth=1.5)+
  geom_vline(xintercept=beta_trait2_dist)+
  geom_vline(xintercept=mean(srmpars$beta_trait2_dist), color="blue", linewidth=1.5, linetype="dashed")+
  xlab("Estimated Distance Effect, Trait 2")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)

latfacresults <- list()
for(i in 1:length(md)){
  tmpmod <- md[[i]][[1]]
  tmpdraws <- extract(tmpmod, pars="latfacterm")
  latfacresults[[i]] <- data.frame(lfmean=apply(tmpdraws$latfacterm, c(2), mean),
                                   lflow=apply(tmpdraws$latfacterm, c(2), quantile, 0.025),
                                   lfhigh=apply(tmpdraws$latfacterm, c(2), quantile, 0.975))
}
for(i in 1:length(md)){
  dtmp <- cbind(md[[i]][[2]], latfacresults[[i]])
  dtmp$rep <- i
  if(i==1){
    d <- dtmp
  }
  if(i>1){
    d <- rbind(d, dtmp)
  }
}

latfacterm_plot <- ggplot(d, aes(x=latfacterm, xend=latfacterm, y=lflow, yend=lfhigh))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  geom_segment(color="blue", alpha=0.1)+
  xlab("True Latent Factor Term")+
  ylab("Estimated Latent Factor Term")+
  theme_bw()+
  theme(aspect.ratio = 1)

library(patchwork)
intercept_plot+sigma_plot+sigma_olre_plot+beta_trait1_sp_plot+beta_trait1_dist_plot+
  beta_trait2_sp_plot+beta_trait2_dist_plot+latfacterm_plot
ggsave("../output_files/sbc_ame_fnch1_plot.png", height=10, width=10)

