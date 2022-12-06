#summarize data for "true trees" to plot ROC or PRC for consensus sequences
#to analyze infector probabilities
#This script will merge all replicates per combination of parameter values
#per depth of sampling

library(stringr)
library(caret)
library(HIVepisimAnalysis)
library(DescTools)
library(dplyr)

rm(merge_label_data)
# this function will merge by perc the W (infector probability) calculated to estimate
# ROC curves
# first merge all data
#then estimate FPR and TPR per perc and param value
merge_label_data <- function(list_dirs, code, sampler){

  #common_dir <- "output/vts/W"

  #empty initial dataframe
  all_data <- data.frame()

  for(i in 1:length(list_dirs)){

    #get combination of parameters, replicate and percentage value
    values <- str_split(list_dirs[i], "/")
    params <- values[[1]][10]
    rep <- values[[1]][11]

    if(sampler == "2"){
      perc <- values[[1]][13]
    } else {
      perc <- values[[1]][12]
    }



    print(i)
    print(paste(params, rep, perc, sep = "."))

    #read simulated data
    #filename <- list.files(paste(list_dirs[i], common_dir, sep = "/"),
    #                       pattern = "W_labels_bst", full.names = TRUE)
    filename <- list.files(list_dirs[i], pattern = "W_labels",
                           full.names = TRUE)
    if(length(filename) >= 2){
      filename <- filename[1]
    }

    #this is to test whether filename exists as a file or is missing
    if(length(filename) != 0){
      if (file.exists(filename) == TRUE){

        #load tm and W to estimate whether a pair is a true pair
        #independet of who infected whom

        data_info <- paste(list_dirs[i], "output/vts/W", sep = "/")
        data_info <- list.files(data_info, pattern = "RData",
                                full.names = TRUE)

        load(data_info)



        sim_data <- read.csv(filename)
        sim_data$labels <- as.character(sim_data$labels)
        sim_data$labels <- as.factor(sim_data$labels)
        sim_data["sampler"] <- sampler
        sim_data["param"] <- str_split(params, "_")[[1]][2]
        sim_data["rep"] <- str_split(rep, "_")[[1]][2]
        sim_data["perc"] <- str_split(perc, "_")[[1]][2]
        sim_data["code"] <- code

        sim_data["group"] <- rep(1:nrow(sim_data))

        sim_data["linked"] <- ifelse(sim_data$infectorProbability >= 0.90, "yes", "no")

        trans_pairs <- sim_data %>% group_by(group) %>%
          group_modify(~check_true_pair(.x, tm))

        all_trans <- trans_pairs %>% group_by(group) %>%
          group_modify(~check_true_transmissions(.x, tm))


        #all_data <- rbind(all_data, sim_data)
        all_data <- rbind(all_data, all_trans)
      }
    }

  }

  return(all_data)

}



# TRUE TREES sampler 1----

#list directories  for sampler 1
#get replicates per parameters value

#250 migrants
rm(replicates)
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_250migrants",
                  full.names = TRUE),
                  full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

#all_data to build ROC or PCR curves
all_data_s1_250mig <- merge_label_data(replicates,
                                       code = "True trees",
                                       sampler = "1")
all_data_s1_250mig["param_perc"] <- paste(all_data_s1_250mig$param,
                                          all_data_s1_250mig$perc,
                                          sep = "_")

saveRDS(all_data_s1_250mig, "all_data_s1_250mig_threshold0.90.RDS")




#500 migrants
rm(replicates)
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "perc")


all_data_s1_500mig <- merge_label_data(replicates,
                                       code = "True trees",
                                       sampler = "1")
all_data_s1_500mig["param_perc"] <- paste(all_data_s1_500mig$param,
                                          all_data_s1_500mig$perc,
                                          sep = "-")

saveRDS(all_data_s1_500mig, "all_data_s1_500mig_threshold0.90.RDS")





#750 migrants
rm(replicates)
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_750migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

all_data_s1_750mig <- merge_label_data(replicates,
                                       code = "True trees",
                                       sampler = "1")
all_data_s1_750mig["param_perc"] <- paste(all_data_s1_750mig$param,
                                          all_data_s1_750mig$perc,
                                          sep = "_")

saveRDS(all_data_s1_750mig, "all_data_s1_750mig_threshold0.90.RDS")




# TRUE TREES sampler 2----

#list directories  for sampler 1
#get replicates per parameters value

#250 migrants
rm(replicates)
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_250migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")

all_data_s2_250mig <- merge_label_data(replicates,
                                      code = "True trees",
                                      sampler = "2")
all_data_s2_250mig["param_perc"] <- paste(all_data_s2_250mig$param,
                                          all_data_s2_250mig$perc,
                                          sep = "_")

saveRDS(all_data_s2_250mig, "all_data_s2_250mig_threshold0.90.RDS")





#500 migrants
rm(replicates)
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")

all_data_s2_500mig <- merge_label_data(replicates,
                                       code = "True trees",
                                       sampler = "2")
all_data_s2_500mig["param_perc"] <- paste(all_data_s2_500mig$param,
                                          all_data_s2_500mig$perc,
                                          sep = "_")

saveRDS(all_data_s2_500mig, "all_data_s2_500mig_threshold0.90.RDS")





#750 migrants
rm(replicates)
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_750migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")


all_data_s2_750mig <- merge_label_data(replicates,
                                       code = "True trees",
                                       sampler = "2")
all_data_s2_750mig["param_perc"] <- paste(all_data_s2_750mig$param,
                                          all_data_s2_750mig$perc,
                                          sep = "_")

saveRDS(all_data_s2_750mig, "all_data_s2_750mig_threshold0.90.RDS")




