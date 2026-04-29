options(scipen=999)
library(compnet)
library(DHARMa)

sim_com_output <- readRDS("sim_com_output.RDS")


output <- data.frame(model=c(),
                     diagnostics=c(),
                     runtime=c())

# Set search parameters and initial dummy values before we search for a data set
# that has reasonable replication of both species and multi-species sites, as well
# as at least some absences.
rich <- 0
coocc_sites <- 0
mat <- 1
## Random data set:
set.seed(181855) # # #
while(rich<10 | coocc_sites<10 | mean(mat)==1){
  g <- sample(c(1:100), size=1)

  mat <- sim_com_output[[g]][[1]]
  mat[mat>0] <- 1
  sums <- rowSums(mat)
  coocc_sites <- length(subset(sums, sums>1))
  rich <- colSums(mat)
  rich[rich>1] <- 1
  rich <- sum(rich)
}
sp_df <- sim_com_output[[g]][[2]]
colnames(mat) <- sp_df$sp
mat <- mat[, colSums(mat)>0]
spp <- colnames(mat)
sp_df <- sp_df[sp_df$sp%in%spp, ]
rownames(sp_df) <- sp_df$sp


## FNCH models without OLRE
# Rank 0:
mod0 <- buildcompnet(presabs = mat,
                      spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                      rank=0,
                      prior_betas_scale=1,
                      iter=2000,
                     family='fnchypg',
                     adapt_delta=0.9,
                      olre=FALSE)
ppred0 <- postpredsamp(mod0)
fpr0 <- apply(ppred0, 1, mean)
dharma0 <- createDHARMa(simulatedResponse = ppred0,
                        observedResponse = mod0$d$both,
                        fittedPredictedResponse = fpr0,
                        integerResponse = TRUE)
testDispersion(dharma0)
testQuantiles(dharma0)
testUniformity(dharma0)
testZeroInflation(dharma0)


# Rank 1
mod1 <- buildcompnet(presabs = mat,
                      spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                      rank=1,
                      prior_betas_scale=1,
                      iter=5000,
                     family='fnchypg',
                     adapt_delta=0.95,
                      olre=FALSE)
ppred1 <- postpredsamp(mod1) # Breaks BiasedUrn; degenerate

output[1, "model"] <- "fnch"
output[1, "diagnostics"] <- "poor"
output[1, "runtime"] <- sum(rstan::get_elapsed_time(mod1$stanmod))


## FNCH models with OLRE
# Rank 0:
mod0 <- buildcompnet(presabs = mat,
                      spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                      rank=0,
                      prior_betas_scale=1,
                      iter=2000,
                     family='fnchypg',
                     adapt_delta=0.9)
ppred0 <- postpredsamp(mod0)
fpr0 <- apply(ppred0, 1, mean)
dharma0 <- createDHARMa(simulatedResponse = ppred0,
                        observedResponse = mod0$d$both,
                        fittedPredictedResponse = fpr0,
                        integerResponse = TRUE)
testDispersion(dharma0)
testQuantiles(dharma0)
testUniformity(dharma0)
testZeroInflation(dharma0)


# Rank 1
mod1 <- buildcompnet(presabs = mat,
                      spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                      rank=1,
                      prior_betas_scale=1,
                      iter=3000,
                     family='fnchypg',
                     adapt_delta=0.95)
ppred1 <- postpredsamp(mod1)
fpr1 <- apply(ppred1, 1, mean)
dharma1 <- createDHARMa(simulatedResponse = ppred1,
                        observedResponse = mod1$d$both,
                        fittedPredictedResponse = fpr1,
                        integerResponse = TRUE)
testDispersion(dharma1) # Breaks BiasedUrn; degenerate

output[2, "model"] <- "fnch_olre"
output[2, "diagnostics"] <- "poor"
output[2, "runtime"] <- sum(rstan::get_elapsed_time(mod1$stanmod))


## Binomial models
# Rank 0:
mod0 <- buildcompnet(presabs = mat,
                     spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                     rank=0,
                     prior_betas_scale=1,
                     iter=2000,
                     adapt_delta=0.9,
                     family="binomial")
ppred0 <- postpredsamp(mod0)
fpr0 <- apply(ppred0, 1, mean)
dharma0 <- createDHARMa(simulatedResponse = ppred0,
                        observedResponse = mod0$d$both,
                        fittedPredictedResponse = fpr0,
                        integerResponse = TRUE)
testDispersion(dharma0)
testQuantiles(dharma0)
testUniformity(dharma0)
testZeroInflation(dharma0)


mod1 <- buildcompnet(presabs = mat,
                     spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                     rank=1,
                     prior_betas_scale=1,
                     iter=3000,
                     adapt_delta=0.9,
                     family="binomial")
ppred1 <- postpredsamp(mod1)
fpr1 <- apply(ppred1, 1, mean)
dharma1 <- createDHARMa(simulatedResponse = ppred1,
                        observedResponse = mod1$d$both,
                        fittedPredictedResponse = fpr1,
                        integerResponse = TRUE)
testDispersion(dharma1)
testQuantiles(dharma1)
testUniformity(dharma1)
testZeroInflation(dharma1)
gofstats(mod1)


mod2 <- buildcompnet(presabs = mat,
                     spvars_dist_int = sp_df[c("ndtrait", "domtrait")],
                     rank=2,
                     prior_betas_scale=1,
                     iter=3000,
                     adapt_delta=0.9,
                     family="binomial")
ppred2 <- postpredsamp(mod2)
fpr2 <- apply(ppred2, 1, mean)
dharma2 <- createDHARMa(simulatedResponse = ppred2,
                        observedResponse = mod2$d$both,
                        fittedPredictedResponse = fpr2,
                        integerResponse = TRUE)
testDispersion(dharma2)
testQuantiles(dharma2)
testUniformity(dharma2)
testZeroInflation(dharma2)
gofstats(mod2)


output[3, "model"] <- "binomial"
output[3, "diagnostics"] <- "poor"
output[3, "runtime"] <- sum(rstan::get_elapsed_time(mod2$stanmod))

saveRDS(output, "modtest_summary_corr18.RDS") # # #
