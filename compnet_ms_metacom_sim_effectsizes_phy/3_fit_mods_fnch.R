options(scipen=999)
library(foreach)
library(doParallel)
registerDoParallel(cores=10)

# Set the number of species
species <- readRDS("species.RDS")

sim_com_output <- readRDS("sim_com_output.RDS")

set.seed(8329)
compnet_res <- foreach(g=1:100) %dopar% {
  gc()
  library(compnet)
  mat <- sim_com_output[[g]][[1]]
  mat[mat>0] <- 1
  sums <- rowSums(mat)
  coocc_sites <- length(subset(sums, sums>1))
  rich <- colSums(mat)
  rich[rich>1] <- 1
  rich <- sum(rich)
  
  if(coocc_sites==0 | rich<3){
    result <- "insufficient"
  }
  
  if(coocc_sites>0 & rich>2){
    sp_df <- sim_com_output[[g]][[2]]
    tree <- sim_com_output[[g]][[3]]
    colnames(mat) <- sp_df$sp
    mat <- mat[, colSums(mat)>0]
    tree <- ape::keep.tip(tree, tip=colnames(mat))
    
    phydistmat <- cophenetic(tree)
    phydist <- data.frame(
      spAid=rownames(phydistmat)[row(phydistmat)[upper.tri(phydistmat)]],
      spBid=colnames(phydistmat)[col(phydistmat)[upper.tri(phydistmat)]],
      phy_dist=phydistmat[upper.tri(phydistmat)])
    
    spp <- colnames(mat)
    
    mod <- buildcompnet(presabs = mat,
                        pairvars = phydist,
                        rank=2,
                        iter=5000,
                        family='fnchypg',
                        prior_betas_scale=1,
                        adapt_delta=0.99)
    result <- mod$stanmod_samp$beta_dy[,1]
  }
  result
}
saveRDS(compnet_res, "compnet_results_fnch.RDS")
