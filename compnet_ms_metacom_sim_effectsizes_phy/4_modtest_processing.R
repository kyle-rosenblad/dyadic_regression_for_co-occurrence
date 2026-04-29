options(scipen=999)

res <- data.frame(dataset=c(),
                  model=c(),
                  performance=c(),
                  runtime=c())
for(i in 1:10){
  temp <- readRDS(paste("modtest_summary_corr", i, ".RDS", sep=""))
  temp$dataset <- i
  res <- rbind(res, temp)
}

saveRDS(res, "../compnet_ms_metacom_sim_summary/modtest_full_summary_corr.RDS")
