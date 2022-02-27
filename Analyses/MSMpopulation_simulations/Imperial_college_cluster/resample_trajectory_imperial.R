#get all log_weights and resample best trajectory
library(stringr)
library(HIVepisimAnalysis)

list_dirs <- dir("/rds/general/user/fferre15/ephemeral", full.names = T,
                 pattern = "large_msm_pop_sim")

list_dirs <- paste(list_dirs, "rep_1", sep = "/")


log_weights_diag <- data.frame()
log_weights_incid <- data.frame()

for(i in 1:length(list_dirs)){

  log_weight_incid_file <- paste(list_dirs[i], "log_weights_incid_df.RDS", sep = "/")
  log_weight_diag_file <- paste(list_dirs[i], "log_weights_dx_df.RDS", sep = "/")


  param <- str_split(list_dirs[i], "/")[[1]][7]
  replicate <- str_split(list_dirs[i], "/")[[1]][8]

  #check whether log_weights for incidence exists
  if(file.exists(log_weight_incid_file)){
    results_incid <- readRDS(log_weight_incid_file)
    results_incid["param"] <- param
    results_incid["replicate"] <- replicate
    results_incid["rep_param"] <- paste(results_incid$replicate, results_incid$param, sep = "_")

    log_weights_incid <- rbind(log_weights_incid, results_incid)
  }

  #check whether log_weights for diagnosis exists
  if(file.exists(log_weight_diag_file)){
    results_dx <- readRDS(log_weight_diag_file)
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


