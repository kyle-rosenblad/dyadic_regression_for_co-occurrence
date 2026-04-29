library(stringr)
library(compnet)
library(ggplot2)
library(patchwork)

results <- data.frame("Gape_Dist_Effect"=c(),
                      "Body_Length_Dist_Effect"=c(),
                      drainage=c())
drainages <- c("Horse", "Powder", "Laramie", "Cheyenne", "South Platte", "Kansas")
for(j in 1:length(drainages)){
  # Load Kirk et al. (2022) trait data and begin formatting for compnet
  traits <- read.csv("Master_trait_analyses.csv")
  rownames(traits) <- traits$Species
  traits[c("Species", "Code", "Number", "Family", "Genus.species")] <- NULL
  
  # Load Kirk et al. (2022) occurrence data, filter to just Horse Creek drainage,
  # and begin formatting for compnet
  pamat <- read.csv("Functional_Dispersion_Analysis_Data.csv")
  sort(table(pamat$Drainage))
  pamat <- subset(pamat, Drainage==drainages[j])
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
  
  # check correlation
  cor(traits[c("Gape", "Max.TL")])
  
  set.seed(90210)
  mod <- buildcompnet(presabs=pamat2,
                      spvars_dist_int=traits[c("Gape", "Max.TL")],
                      prior_betas_scale=1,
                      iter=3000,
                      family="binomial")
  results2 <- as.data.frame(mod$stanmod_samp$beta_dy)
  names(results2) <- c("Gape_Dist_Effect", "Body_Length_Dist_Effect")
  results2$drainage <- drainages[j]
  results <- rbind(results, results2)
}

p1 <- ggplot()+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_density(data=subset(results, drainage!="Horse"), aes(x=Gape_Dist_Effect, group=drainage), color="gray80")+
  geom_density(data=subset(results, drainage=="Horse"), aes(x=Gape_Dist_Effect, group=drainage), color="red", linewidth=1)+
  xlab("Standardized Gape Width Distance Effect")+
  ylab("Density")+
  theme_bw()

p2 <- ggplot()+
  geom_vline(xintercept=0)+
  geom_hline(yintercept=0)+
  geom_density(data=subset(results, drainage!="Horse"), aes(x=Body_Length_Dist_Effect, group=drainage), color="gray80")+
  geom_density(data=subset(results, drainage=="Horse"), aes(x=Body_Length_Dist_Effect, group=drainage), color="red", linewidth=1)+
  xlab("Standardized Body Length Distance Effect")+
  theme_bw()+
  theme(axis.title.y=element_blank())
  

p1+p2
ggsave("../output_files/compnet_kirk_geb_alternate_drainages.png", width=7.2, height=4.8)
