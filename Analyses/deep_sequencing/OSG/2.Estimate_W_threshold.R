#Script to estimate the best W threshold to be use in subsequent analysis


targz_files <- list.files(list.files(list.files(list.files("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/deepseq",
                         full.names = TRUE, pattern = "best_trajectories"),
                         full.names = TRUE),
                         full.names = TRUE),
                         full.names = TRUE, pattern = "processing_network_results.tar.gz")

real_trans_W_tm <- data.frame()

for (i in 1:length(targz_files)){

  untar(targz_files[i])
  load("output_deepseq/vts/merged_trees_sampling_migrant_years_1_simple_0.9.RData")

  if(length(real_trans$donor_ID) != 0 ){
    params <- str_split(targz_files[i], pattern = "/")

    real_trans["params"] <- params[[1]][11]
    real_trans["rep"] <- params[[1]][12]
    real_trans["migrants"] <- str_split(params[[1]][10], pattern = "_")[[1]][3]



    real_trans_W_tm <- rbind(real_trans_W_tm, real_trans)
  }




}

real_trans_W_tm["params_mig"] <- paste(real_trans_W_tm$params, real_trans_W_tm$migrants, sep = "_")

library(ggplot2)

ggplot(real_trans_W_tm, aes(x=migrants, y=infectorProbability, fill=migrants)) +
  geom_boxplot()

quartz()
ggplot(real_trans_W_tm, aes(x=params_mig, y=infectorProbability, fill=migrants)) +
  geom_boxplot()

mig250 <- subset(real_trans_W_tm, migrants == "250migrants")
mig500 <- subset(real_trans_W_tm, migrants == "500migrants")
mig750 <- subset(real_trans_W_tm, migrants == "750migrants")

mig250_1067 <- subset(real_trans_W_tm, migrants == "250migrants" & params == "params_1067")
mig500_1067 <- subset(real_trans_W_tm, migrants == "500migrants" & params == "params_1067")
mig750_1067 <- subset(real_trans_W_tm, migrants == "750migrants" & params == "params_1067")

mig250_2348 <- subset(real_trans_W_tm, migrants == "250migrants" & params == "params_2348")
mig500_2348 <- subset(real_trans_W_tm, migrants == "500migrants" & params == "params_2348")
mig750_2348 <- subset(real_trans_W_tm, migrants == "750migrants" & params == "params_2348")

ggplot(mig500_1067, aes(x=params_mig, y=infectorProbability, fill=migrants)) +
  geom_boxplot()
