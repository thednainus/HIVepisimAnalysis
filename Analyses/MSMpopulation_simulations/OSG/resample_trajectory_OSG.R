#get all log_weights and resample best trajectory
library(stringr)
library(HIVepisimAnalysis)

#source observed data
#source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))

#list_dirs <- dir("/rds/general/user/fferre15/ephemeral/results_small_pop", full.names = T)

#list_dirs <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/results_params_24Jan22", full.names = T)
list_dirs <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results", full.names = T)
list_dirs <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results_art_2004", full.names = T)
list_dirs <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results_narrow_parameters", full.names = T)


list_dirs <- paste(list_dirs, "rep_1", sep = "/")


log_weights_diag <- data.frame()
log_weights_incid <- data.frame()

for(i in 1:length(list_dirs)){

  log_weight_file <- paste(list_dirs[i], "log_weights.tar.gz", sep = "/")
  untar(log_weight_file)

  #teste <- readRDS("log_weights/log_weights_incid_df.RDS")

  param <- str_split(list_dirs[i], "/")[[1]][8]
  replicate <- str_split(list_dirs[i], "/")[[1]][9]

  #check whether log_weights for incidence exists
  if(file.exists("log_weights/log_weights_incid_df.RDS")){
    results_incid <- readRDS("log_weights/log_weights_incid_df.RDS")
    results_incid["param"] <- param
    results_incid["replicate"] <- replicate
    results_incid["rep_param"] <- paste(results_incid$replicate, results_incid$param, sep = "_")

    log_weights_incid <- rbind(log_weights_incid, results_incid)
  }

  #check whether log_weights for diagnosis exists
  if(file.exists("log_weights/log_weights_dx_df.RDS")){
    results_dx <- readRDS("log_weights/log_weights_dx_df.RDS")
    results_dx["param"] <- param
    results_dx["replicate"] <- replicate
    results_dx["rep_param"] <- paste(results_dx$replicate, results_dx$param, sep = "_")

    log_weights_diag <- rbind(log_weights_diag, results_dx)
  }

}

#resample for incidence
w_incid =   exp( log_weights_incid[,2]  - max( log_weights_incid[,2] ) )
resample_incid <- sample( as.character( log_weights_incid[, 'rep_param']  ),
                          prob = w_incid, replace = TRUE )

#resample for diagnosis
w_dx =   exp( log_weights_diag[,2]  - max( log_weights_diag[,2] ) )
resample_dx <- sample( as.character( log_weights_diag[, 'rep_param']  ),
                       prob = w_dx, replace = TRUE )

saveRDS(unique(resample_incid), "resample_incid.RDS")
saveRDS(log_weights_incid, "log_weights_incid_all_df.RDS")

saveRDS(unique(resample_dx), "resample_diag.RDS")
saveRDS(log_weights_diag, "log_weights_diag_all_df.RDS")








#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

sim_df <- as.data.frame(sim)
#sim_df$sim <- "1"
sim_df$type <- paste("param", as.character(line_number), sep = "_")

pop1_data <- sim_df[c("type", "time", "newDx_pop1", "incid.pop1", "num.pop1", "i.num.pop1")]
pop1_data_freq <- data.frame(type = pop1_data["type"],
                             time = pop1_data["time"],
                             newDx_pop1 = pop1_data["newDx_pop1"],
                             incidence_pop1 = pop1_data["incid.pop1"])


pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

pop1_data <- pop1_data_freq
pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
pop1_data$year <- as.factor(pop1_data$year)
pop1_data$type <- as.factor(pop1_data$type)
#pop1_data$sim <- as.factor(pop1_data$sim)
pop1_data <- pop1_data[c(6, 1:4)]

#aggregate by new HIV diagnosis and by incidence for population 1
newDx_pop1_agg <- aggregate(newDx_pop1  ~ year + type, data = pop1_data, FUN=sum)
incid_pop1_agg <- aggregate(incid.pop1  ~ year + type, data = pop1_data, FUN=sum)


newDx_pop1_agg["rep_param"] <- newDx_pop1_agg$type
incid_pop1_agg["rep_param"] <- incid_pop1_agg$type

library(ggplot2)
ggplot(newDx_pop1_agg[c(6:20),], aes(x=year)) +
  geom_line(aes(y = newDx_pop1)) +
  theme_bw()

ggplot(incid_pop1_agg[c(6:20),], aes(x=year)) +
  geom_line(aes(y = incid.pop1)) +
  theme_bw()

#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence_model <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                       package = "HIVepisimAnalysis"))



function ( diag_obs, diag_sim )
{
  #browser()
  diag_sim <- diag_sim$newDx_pop1
  #subset to remove the NA in observed data (data not available from 1981 and 2021)
  #diag_sim <- diag_sim[3:41]
  #diag_obs <- diag_obs[3:41]

  diag_sim <- diag_sim[6:41]
  diag_obs <- diag_obs[6:41]

  log_importance_weight <- sum( dpois( diag_sim, lambda = diag_obs , log = TRUE ) )

  return(log_importance_weight)
}

#new HIV diagnosis
log_weights_dx <- compute_log_importance_weight_newDx(incidenceDiag$frequency,
                                                      diag_sim = newDx_pop1_agg)
log_weights_dx_df <- data.frame(rep_param = newDx_pop1_agg$rep_param[1],
                                log_weights = log_weights_dx)

#new HIV incidence
log_weights_incid <- compute_log_importance_weight_incidence(incidence_model$ECDC_incidence.N_Inf_M,
                                                             incid_sim = incid_pop1_agg[1:41,])
log_weights_incid_df <- data.frame(rep_param = incid_pop1_agg$rep_param[1], log_weights = log_weights_incid)


saveRDS(log_weights_dx_df, "log_weights_dx_df.RDS")
saveRDS(newDx_pop1_agg, "summary_newDx_pop1.RDS")


saveRDS(log_weights_incid_df, "log_weights_incid_df.RDS")
saveRDS(incid_pop1_agg, "summary_incidence_pop1.RDS")


