library(ggplot2)
library(compnet)

species <- 15
sites <- 1005

set.seed(392)
bl <- rep(c(-2:2), times=species/5)
gw <- bl+rnorm(n=species, sd=0.5)

traits <- data.frame(bl, gw)
traits$sp <- c(1:nrow(traits))
rownames(traits) <- traits$sp

env <- rep_len(gw, length.out=sites)

spsite <- as.data.frame(merge(c(1:species), c(1:sites)))
names(spsite) <- c("sp", "site")
spsite <- spsite[sample(nrow(spsite)),]

presabs <- matrix(0, nrow=sites, ncol=species)

for(k in 1:nrow(spsite)){
  i <- spsite[k, "sp"]
  j <- spsite[k, "site"]
  
  occprob <- 1
  
  envdist <- abs(env[j]-traits[i, "bl"])
  envfac <- exp(-envdist)
  
  occprob <- occprob*envfac
  
  focalgw <- gw[i]
  otherspp <- presabs[j,]
  othergw <- gw[which(otherspp==1)]
  othergw <- othergw[othergw!=focalgw]
  if(length(othergw)>0){
    compdists <- 0.5*abs(othergw-focalgw)
    compfacs <- exp(-compdists)
    compfac <- sum(compfacs)
    occprob <- occprob-compfac
    if(occprob<0){
      occprob <- 0
    }
  }
  presabs[j,i] <- rbinom(n=1, prob=occprob, size=1)
}

colnames(presabs) <- c(1:ncol(presabs))
rownames(presabs) <- c(1:nrow(presabs))

# quick dummy compnet model to get needed data structure
mod_dist <- buildcompnet(presabs=presabs,
                         spvars_dist_int=traits[c("bl", "gw")],
                         iter=2,
                         warmup=1)
sppairs <- as.data.frame(mod_dist$d)

plot1 <- ggplot(sppairs, aes(x=gw_dist, y=both/either, color=bl_dist, succ=both, fail=either-both))+
  geom_point()+
  geom_smooth(
    color="black",
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(succ, fail) ~ x)+
  ggtitle("No Confounding Control")+
  xlab("Gape Width Distance")+
  ylab("P(Co-occurrence)")+
  scale_color_viridis_c(name="Body\nLength\nDistance",
                        option="plasma")+
  theme_bw()+
  theme(aspect.ratio=1,
        legend.position="none")

plot2 <- ggplot(sppairs, aes(x=gw_dist, y=both/either, succ=both, fail=either-both, group=bl_dist, color=bl_dist, fill=bl_dist))+
  geom_point()+
  stat_smooth(
    method="glm",
    method.args=list(family="binomial"),
    formula = cbind(succ, fail) ~ x)+
  scale_color_viridis_c(name="Body\nLength\nDistance",
                        option="plasma")+
  scale_fill_viridis_c(name="Body\nLength\nDistance",
                       option="plasma")+
  ggtitle("With Confounding Control")+
  xlab("Gape Width Distance")+
  ylab("P(Co-occurrence)")+
  theme_bw()+
  theme(aspect.ratio=1,
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
library(patchwork)
plot1+plot2
ggsave("../output_files/fig1_scatterplots.pdf", height=3, width=6)
