options(scipen=999)
library(ggplot2)
library(ggthemes)
library(patchwork)

res <- data.frame(prop_pos_bin=c())

r <- readRDS("compnet_results_bin.RDS")
for(i in 1:length(r)){
  if(class(r[[i]])!="character"){
    res[i, "prop_pos_bin"] <- mean(r[[i]]>0)
  }
  if(class(r[[i]])=="character"){
    res[i, "prop_pos_bin"] <- NA
  }
}
table(res$prop_pos_bin>0.975)
table(res$prop_pos_bin<0.025)



rk0 <- readRDS("compnet_results_bin_k0.RDS")
for(i in 1:length(rk0)){
  if(class(rk0[[i]])!="character"){
    res[i, "prop_pos_bin_k0"] <- mean(rk0[[i]]>0)
  }
  if(class(rk0[[i]])=="character"){
    res[i, "prop_pos_bin_k0"] <- NA
  }
}
table(res$prop_pos_bin_k0>0.975)
table(res$prop_pos_bin_k0<0.025)




r2 <- readRDS("compnet_results_fnch.RDS")
for(i in 1:length(r2)){
  if(class(r2[[i]])!="character"){
    res[i, "prop_pos_fnch"] <- mean(r2[[i]]>0)
  }
  if(class(r2[[i]])=="character"){
    res[i, "prop_pos_fnch"] <- NA
  }
}
table(res$prop_pos_fnch>0.975)
table(res$prop_pos_fnch<0.025)





r2k0 <- readRDS("compnet_results_fnch_k0.RDS")
for(i in 1:length(r2k0)){
  if(class(r2k0[[i]])!="character"){
    res[i, "prop_pos_fnch_k0"] <- mean(r2k0[[i]]>0)
  }
  if(class(r2[[i]])=="character"){
    res[i, "prop_pos_fnch_k0"] <- NA
  }
}
table(res$prop_pos_fnch_k0>0.975)
table(res$prop_pos_fnch_k0<0.025)





