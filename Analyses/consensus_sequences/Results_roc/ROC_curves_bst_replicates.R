#ROC curves for analyzing infector probabilities
#This script will merge all replicates per combination of parameter values
#per depth of sampling

library(stringr)
library(caret)
library(HIVepisimAnalysis)
library(DescTools)


# this function will merge by perc the W (infector probability) calculated to estimate
# ROC curves
# first merge all data
#then estimate FPR and TPR per perc and param value
merge_label_data <- function(list_dirs, code, sampler){

  common_dir <- "output/vts/W"

  #empty initial dataframe
  all_data <- data.frame()

  for(i in 1:length(list_dirs)){

    #get combination of parameters, replicate and percentage value
    values <- str_split(list_dirs[i], "/")
    params <- values[[1]][6]
    rep <- values[[1]][7]
    perc <- values[[1]][8]

    print(i)
    print(paste(params, rep, perc, sep = "."))

    #read simulated data
    filename <- list.files(paste(list_dirs[i], common_dir, sep = "/"),
                           pattern = "W_labels_bst", full.names = TRUE)


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

replicates <- dir(dir(dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_250migrants",
                  full.names = TRUE),
                  full.names = TRUE),
                  full.names = TRUE, pattern = "perc")

all_roc_data_s1_50mig <- merge_label_data(replicates,
                                          code = "True trees",
                                          sampler = "1")
all_roc_data_s1_50mig["param_perc"] <- paste(all_roc_data_s1_50mig$param,
                                             all_roc_data_s1_50mig$perc,
                                             sep = "_")

all_roc_data_s1_50mig <- readRDS("all_roc_label_data_s1_50mig.RDS")

get_roc_data <- function(df_data){

  thresholds <- seq(from = 0, to = 1, by = 0.01)
  rates <- lapply(thresholds, get_rates2, df_data)
  rates <- do.call(rbind, rates)

  return(rates)

}


roc_data <- all_roc_data_s1_50mig %>%
  group_by(param_perc) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], sep = "_"))
    #print(unique(.x$rep))
    get_roc_data(.x)

  })

teste_data1 <- subset(all_roc_data_s1_50mig,
                     param_perc == "1067_0.05")

teste_roc1 <- get_roc_data(teste_data1)

teste_data2 <- subset(all_roc_data_s1_50mig,
                      param_perc == "2348_0.05")
teste_roc2 <- get_roc_data(teste_data1)

teste_roc <- do.call(rbind,roc_data)

all_true_roc_s1 <- merge_roc_data(replicates,
                                  code = "True trees",
                                  sampler = 1)

saveRDS(all_true_roc_s1, "all_true_roc_50migrants_sampler1_bst.RDS")






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



