library(stringr)
library(compnet)
library(ggplot2)
library(patchwork)


# Load Kirk et al. (2022) trait data and begin formatting for compnet
traits <- read.csv("Master_trait_analyses.csv")
rownames(traits) <- traits$Species
traits[c("Species", "Code", "Number", "Family", "Genus.species")] <- NULL

# Load Kirk et al. (2022) occurrence data, filter to just Horse Creek drainage,
# and begin formatting for compnet
pamat <- read.csv("Functional_Dispersion_Analysis_Data.csv")
pamat <- subset(pamat, Drainage=="Horse")
uniquesites <- unique(pamat$Repeat)
# Create a new presence-absence matrix for unique spatial locations
for(i in 1:length(uniquesites)){
  tmp <- subset(pamat, Repeat==uniquesites[i])
  tmp <- tmp[1,]
  if(i==1){
    pamat2 <- tmp
  }
  if(i>1){
    pamat2 <- rbind(pamat2, tmp)
  }
}
rownames(pamat2) <- pamat2$Repeat # make rownames unique site names
pamat2 <- pamat2[,16:ncol(pamat2)]
pamat2 <- as.matrix(pamat2)
pamat2 <- pamat2[,colSums(pamat2)>0] # drop species absent from this drainage
colnames(pamat2) <- str_replace_all(colnames(pamat2), "\\.", " ") # fix species names to match trait data

traits <- traits[rownames(traits)%in%colnames(pamat2),] # filter irrelevant species from trait data
pamat2 <- pamat2[,colnames(pamat2)%in%rownames(traits)] # drop species with no trait data

# Show correlation
cor(traits[c("Gape", "Max.TL")])

set.seed(90210)
mod <- buildcompnet(presabs=pamat2,
                    spvars_dist_int=traits[c("Gape", "Max.TL")],
                    prior_betas_scale=1,
                    iter=3000,
                    family="binomial")
mod2 <- buildcompnet(presabs=pamat2,
                     spvars_dist_int=traits[c("Gape")],
                     prior_betas_scale=1,
                     iter=3000,
                     family="binomial")
mod3 <- buildcompnet(presabs=pamat2,
                     spvars_dist_int=traits[c("Max.TL")],
                     prior_betas_scale=1,
                     iter=3000,
                     family="binomial")
summarize_compnet(mod)
summarize_compnet(mod2)
summarize_compnet(mod3)

p1 <- ggplot(traits, aes(x=Max.TL, y=Gape))+
  annotate("text", label = "a", x = -Inf, y = Inf, size=10, hjust=-0.5, vjust=1.5)+
  geom_point(size=3)+
  xlab("Body Length (cm)")+
  ylab("Gape Width (mm)")+
  theme_bw()+
  theme(aspect.ratio=1)
p1

gape_effects2 <- data.frame(es=mod$stanmod_samp$beta_dy[,1],
                          Model="Two-Trait")
gape_effects1 <- data.frame(es=mod2$stanmod_samp$beta_dy[,1],
                           Model="One-Trait")
body_effects2 <- data.frame(es=mod$stanmod_samp$beta_dy[,2],
                             Model="Two-Trait")
body_effects1 <- data.frame(es=mod3$stanmod_samp$beta_dy[,1],
                             Model="One-Trait")

gape_effects <- rbind(gape_effects1, gape_effects2)
body_effects <- rbind(body_effects1, body_effects2)


p2 <- ggplot(gape_effects, aes(x=es, group=Model, color=Model))+
  annotate("text", label = "b", x = -Inf, y = Inf, size=10, hjust=-0.5, vjust=1.5)+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_density(size=1)+
  scale_color_manual(values=c("red", "blue"))+
  xlab("Effect of Gape Width Distance on Co-occurrence")+
  ylab("Estimated Posterior Density")+
  theme_bw()+
  theme(aspect.ratio=1)
p2


p3 <- ggplot(body_effects, aes(x=es, group=Model, color=Model))+
  annotate("text", label = "c", x = -Inf, y = Inf, size=10, hjust=-0.5, vjust=1.5)+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_density(size=1)+
  scale_color_manual(values=c("red", "blue"))+
  xlab("Effect of Body Length Distance on Co-occurrence")+
  ylab("Estimated Posterior Density")+
  theme_bw()+
  theme(aspect.ratio=1)
p3




p1+p2+p3+plot_layout(guides="collect")
ggsave("../output_files/fig4.pdf", width=12, height=4)
ggsave("../output_files/fig4.png", width=12, height=4)
