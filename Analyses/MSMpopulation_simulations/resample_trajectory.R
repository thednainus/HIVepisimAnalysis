#get all log_weights and resample best trajectory
library(stringr)
library(HIVepisimAnalysis)

#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))

list_dirs <- dir("/rds/general/user/fferre15/ephemeral/results_small_pop", full.names = T)

log_weights <- data.frame()

for(i in 1:length(list_dirs)){

    for(j in 1:50){
      new_dir=paste(list_dirs[i], "/rep_", j, sep ="")

      param_rep <- str_split(new_dir, "/")

      if(file.exists(paste(new_dir, "log_importance_wghts.RDS", sep = "/"))){
        results <- readRDS(paste(new_dir, "log_importance_wghts.RDS", sep = "/"))
        results["param"] <- param_rep[[1]][8]
        results["param"] <- unlist(lapply(results$param, function(x) str_replace(string = x,
                                                                                 pattern = "small_msm_pop_sim",
                                                                                 replacement = "param")))
        results["replicate"] <- param_rep[[1]][9]
        results["rep_param"] <- paste(results$replicate, results$param, sep = "_")

        log_weights <- rbind(log_weights, results)
      }



    }
}


w =   exp( log_weights[,2]  - max( log_weights[,2] ) )
resample <- sample( as.character( log_weights[, 'rep_param']  ), prob = w, replace = TRUE )

saveRDS(resample, "resample.RDS")
saveRDS(log_weights, "log_weights_all_df.RDS")
