options(scipen=999)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(officer)
library(flextable)
library(brms)

res <- data.frame(dataset=c(),
                  model=c(),
                  performance=c(),
                  runtime=c())
for(i in 1:20){
  temp <- readRDS(paste("modtest_summary_corr", i, ".RDS", sep=""))
  temp$dataset <- i
  res <- rbind(res, temp)
}

res$diagnostics <- factor(res$diagnostics, levels=c("good", "fair", "poor"))
table(res[c("model", "diagnostics")])


runtimemod <- brm(runtime ~ model + (1|dataset),
                  family=Gamma(link="log"),
                  data=subset(res, model!="fnch"),
                  chains=2,
                  cores=2)
summary(runtimemod)
runtimemod_eff <- conditional_effects(runtimemod)$model
runtimemod_eff

runtime_plot <- ggplot(runtimemod_eff)+
  geom_point(aes(x=model, y=estimate__), size=3)+
  geom_segment(aes(x=model, y=lower__, yend=upper__))+
  theme_bw()+
  scale_x_discrete(labels=c("Binomial", "FNCH OLRE"))+
  xlab("Model Version")+
  ylab("Estimated Runtime (s)")+
  theme(aspect.ratio=1)
runtime_plot
ggsave("../output_files/runtime_plot.png", height=3, width=3)
ggsave("../output_files/runtime_plot.pdf", height=3, width=3)
