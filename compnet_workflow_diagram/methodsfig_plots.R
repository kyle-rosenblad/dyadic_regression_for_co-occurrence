library(ggplot2)

set.seed(2390)
dgape <- data.frame(effectsize = rnorm(n=3e4, mean=2),
                trait = "Gape Width")
dbody <- data.frame(effectsize = rnorm(n=3e4, mean=-2),
                trait = "Body Length")
d <- rbind(dgape, dbody)
ggplot(d, aes(x=effectsize, color=trait))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_density(linewidth=1)+
  scale_color_manual(values=c("orange", "purple"), name="Trait")+
  ggtitle("Trait Distance Effects on Cooccurrence")+
  xlab("Posterior Effect Size Estimate")+
  ylab("Density")+
  theme_bw()+
  theme(aspect.ratio=1)
ggsave("../output_figures/MethodsFigPlot.png", height=3, width=4)
