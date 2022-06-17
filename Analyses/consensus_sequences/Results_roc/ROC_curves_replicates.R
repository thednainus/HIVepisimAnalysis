#ROC curves for analyzing infector probabilities
#This script will merge all replicates per combination of parameter values
#per depth of sampling

library(stringr)
library(caret)
library(HIVepisimAnalysis)
library(DescTools)
library(dplyr)


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


    #this is to test whether filename exists as a file or is missing
    if(length(filename) != 0){
      if (file.exists(filename) == TRUE){

        sim_data <- read.csv(filename)
        sim_data$labels <- as.character(sim_data$labels)
        sim_data$labels <- as.factor(sim_data$labels)
        sim_data["sampler"] <- sampler
        sim_data["param"] <- str_split(params, "_")[[1]][2]
        sim_data["rep"] <- str_split(rep, "_")[[1]][2]
        sim_data["perc"] <- str_split(perc, "_")[[1]][2]
        sim_data["code"] <- code


        all_data <- rbind(all_data, sim_data)
      }
    }

  }

  return(all_data)

}



# TRUE TREES sampler 1----

#list directories  for sampler 1
#get replicates per parameters value

#250 migrants
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_250migrants",
                  full.names = TRUE),
                  full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s1_250mig <- merge_label_data(replicates,
                                          code = "True trees",
                                          sampler = "1")
all_roc_data_s1_250mig["param_perc"] <- paste(all_roc_data_s1_250mig$param,
                                             all_roc_data_s1_250mig$perc,
                                             sep = "_")

saveRDS(all_roc_data_s1_250mig, "all_roc_data_s1_250mig.RDS")



#500 migrants
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s1_500mig <- merge_label_data(replicates,
                                           code = "True trees",
                                           sampler = "1")
all_roc_data_s1_500mig["param_perc"] <- paste(all_roc_data_s1_500mig$param,
                                              all_roc_data_s1_500mig$perc,
                                              sep = "_")

saveRDS(all_roc_data_s1_500mig, "all_roc_data_s1_500mig.RDS")




#750 migrants
replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_750migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s1_750mig <- merge_label_data(replicates,
                                           code = "True trees",
                                           sampler = "1")
all_roc_data_s1_750mig["param_perc"] <- paste(all_roc_data_s1_750mig$param,
                                              all_roc_data_s1_750mig$perc,
                                              sep = "_")

saveRDS(all_roc_data_s1_750mig, "all_roc_data_s1_750mig.RDS")



# TRUE TREES sampler 2----

#list directories  for sampler 1
#get replicates per parameters value

#250 migrants
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_250migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s2_250mig <- merge_label_data(replicates,
                                           code = "True trees",
                                           sampler = "2")
all_roc_data_s2_250mig["param_perc"] <- paste(all_roc_data_s2_250mig$param,
                                              all_roc_data_s2_250mig$perc,
                                              sep = "_")

saveRDS(all_roc_data_s2_250mig, "all_roc_data_s2_250mig.RDS")



#500 migrants
rm(replicates)
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s2_500mig <- merge_label_data(replicates,
                                           code = "True trees",
                                           sampler = "2")
all_roc_data_s2_500mig["param_perc"] <- paste(all_roc_data_s2_500mig$param,
                                              all_roc_data_s2_500mig$perc,
                                              sep = "_")

saveRDS(all_roc_data_s2_500mig, "all_roc_data_s2_500mig.RDS")




#750 migrants
rm(replicates)
replicates <- dir(dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_750migrants",
                          full.names = TRUE),
                      full.names = TRUE),
                  full.names = TRUE, pattern = "sampler2"),
                  full.names = TRUE, pattern = "perc")


all_roc_data_s2_750mig <- merge_label_data(replicates,
                                           code = "True trees",
                                           sampler = "2")
all_roc_data_s2_750mig["param_perc"] <- paste(all_roc_data_s2_750mig$param,
                                              all_roc_data_s2_750mig$perc,
                                              sep = "_")

saveRDS(all_roc_data_s2_750mig, "all_roc_data_s2_750mig.RDS")



#read for ROC curves
all_roc_data_s1_250mig <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/all_roc_data_s1_250mig.RDS")
#all_roc_data_s1_500mig <- readRDS("all_roc_data_s1_500mig.RDS")
#all_roc_data_s1_750mig <- readRDS("all_roc_data_s1_750mig.RDS")


get_roc_data <- function(df_data){

  thresholds <- seq(from = 0, to = 1, by = 0.01)
  rates <- lapply(thresholds, get_rates2, df_data)
  rates <- do.call(rbind, rates)
  rates <- distinct(rates)

  return(rates)

}

all_roc_data_s1_250mig["param_perc_rep"] <- paste(all_roc_data_s1_250mig$param,
                                                  all_roc_data_s1_250mig$perc,
                                                  all_roc_data_s1_250mig$rep,
                                                  sep = "_")

teste1.1 <- subset(all_roc_data_s1_250mig, rep == "1" & param == "2348" &
                     perc == "0.05")
teste1.30 <- subset(all_roc_data_s1_250mig, rep == "30" & param == "1067" & perc == "0.05")
teste1.1.roc <- get_roc_data(teste1.1)
teste1.30.roc <- get_roc_data(teste1.30)

roc_data <- all_roc_data_s1_250mig %>%
  group_by(param_perc_rep) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], sep = "_"))
    #print(unique(.x$rep))
    get_roc_data(.x)

  })
teste2=do.call(rbind, roc_data)
teste2.1.roc <- subset(teste2, param == "1067" & rep == "1" & perc == "0.05" )
teste2.30.roc <- subset(teste2, param == "1067" & rep == "30" & perc == "0.05" )

#saveRDS(roc_data, "roc_data_250mig.RDS")
#saveRDS(roc_data, "roc_data_500mig.RDS")
saveRDS(roc_data, "roc_data_750mig.RDS")




#random sample 5,000 points for plotting
random_data <- all_roc_data_s1_50mig %>%
  group_by(param_perc) %>%
  group_map(~ {
    .x[sample(nrow(.x), 5000), ]


  })



auc <- all_true_roc_s1 %>%
  group_by(param_perc) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })

auc2 <- all_true_roc_s1_noNA %>%
  group_by(param_perc) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })


quartz()
ggplot(data = all_true_roc_s1_noNA, aes(x=FPR, y=TPR)) +
  geom_line(size = 1) +
  geom_abline(slope= 1, intercept= 0, linetype = 4) +
  theme_bw(base_size = 20) +
  facet_wrap(~ param_perc, scales = "free", ncol = 4)+
  theme(legend.position="none") +
  ggtitle("ROC curve") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate")



