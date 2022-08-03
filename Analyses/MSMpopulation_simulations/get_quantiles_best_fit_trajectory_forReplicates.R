#plot best fit trajectory
library(EpiModel)
library(HIVepisimAnalysis)
library(lubridate)
library(stringr)


m_and_qt <- function(dataframe){

  all_qt <- aggregate(dataframe[,3] ~ year, data = dataframe,
                      FUN=quantile, probs = c(0.025, 0.5, 0.975))

  all_qt <- data.frame(year = all_qt$year, all_qt$`dataframe[, 3]`)

  return(all_qt)


}

#250 migrants
dir_list <- dir("../Results_paper/best_trajectories_250migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_250migrants/params_2348", full.names = T)

#500 migrants
dir_list <- dir("../Results_paper/best_trajectories_500migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_500migrants/params_2348", full.names = T)


#750 migrants
dir_list <- dir("../Results_paper/best_trajectories_75migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_75migrants/params_2348", full.names = T)




#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

results_incid <- data.frame()
results_diag <- data.frame()

for(i in 1:length(dir_list)){



  #untar results files
  results_tar_files <- list.files(dir_list[i], full.names = T,
                                  pattern = "sim_results.tar.gz")
  untar(results_tar_files)

  filename <- "results/results_sim.RDS"

  sim <- readRDS(filename)

  sim_param <- str_split(results_tar_files, "/")

  sim_df <- as.data.frame(sim)
  sim_df$param <- str_split(sim_param[[1]][7], "/")[[1]][1]

  sim_df$replicate <- sim_param[[1]][8]

  sim_df["rep_param"] <- paste(sim_df$replicate, sim_df$param, sep = "_")

  #diagnosis
  pop1_data_freq <- sim_df[c("rep_param", "time", "newDx_pop1")]
  pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

  pop1_data <- pop1_data_freq
  pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
  pop1_data$year <- as.factor(pop1_data$year)
  pop1_data$type <- as.factor(pop1_data$rep_param)

  pop1_data <- pop1_data[c(5, 1,3)]
  pop1_data_agg <- aggregate(newDx_pop1 ~ year + rep_param, data = pop1_data, FUN=sum)

  #incidence
  pop1_data_incid <- sim_df[c("rep_param", "time", "incid.pop1")]
  pop1_data_incid["time_years"] <- days2years(pop1_data_incid$time, init_date = init_sim_date)

  pop1_data_incid2 <- pop1_data_incid
  pop1_data_incid2["year"] <- unlist(lapply(pop1_data_incid2$time_years,
                                            function(x) strsplit(as.character(x),
                                                                 split = ".", fixed = TRUE)[[1]][1]))
  pop1_data_incid2$year <- as.factor(pop1_data_incid2$year)
  pop1_data_incid2$type <- as.factor(pop1_data_incid2$rep_param)

  pop1_data_incid2 <- pop1_data_incid2[c(5, 1,3)]
  pop1_data_agg_incid <- aggregate(incid.pop1 ~ year + rep_param, data = pop1_data_incid2, FUN=sum)

  results_diag <- rbind(results_diag, pop1_data_agg)
  results_incid <- rbind(results_incid, pop1_data_agg_incid)

}

#get median and upper and lower quantiles
all_diag <- m_and_qt(results_diag)
all_incid <- m_and_qt(results_incid)

saveRDS(all_diag, "all_diag_m_and_q_2348_500migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_500migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_500migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_500migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_meandegree.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_75migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_75migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_25migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_25migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_meandegree.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_75migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_75migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_25migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_25migrants.RDS")
