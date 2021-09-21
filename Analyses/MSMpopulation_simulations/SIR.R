#Sampling importance resampling (SIR) using deviance for poisson distribtution

#deviance for poisson distribution

pop1_data_agg
pop1_data_agg["rep_param"] <- paste(pop1_data_agg$param, "replicate", pop1_data_agg$sim, sep = "_")
pop1_data_agg$rep_param <- as.factor(pop1_data_agg$rep_param)
pop1_data_agg <- pop1_data_agg[,c(1,6,4)]

#sampling importance resampling
test <- ddply(pop1_data_agg, "rep_param", sir,
               diag_obs = incidenceDiag$frequency)


sir <- function(diag_obs, diag_est){
  #browser()
  diag_est <- diag_est$newDx_pop1
  #subset to remove the NA in observed data (data not available fro 1981 and 2021)
  diag_est <- diag_est[3:41]

  diag_obs <- diag_obs[3:41]

  dev_poi <- list()

  for(i in 1:length(diag_obs)){
    #I used the equation for the unit deviance for the poisson distribution
    dev_poi[[i]] <- 2 * (diag_obs[i] * log(diag_obs[i]/diag_est[i]) - diag_obs[i] + diag_est[i])

  }

  dev_poi <- unlist(dev_poi)
  #browser()
  #normalize deviance
  dev_poi <- dev_poi/sum(dev_poi)

  browser()

  if(sum(dev_poi) > 0 & is.na(sum(dev_poi)) == FALSE ){
    new_sampling <- sample(x = diag_est,
                           size = 39,
                           replace = TRUE,
                           prob = dev_poi)
  }else{
    new_sampling <- rep(NA, 39)
  }

  return(new_sampling)
}
