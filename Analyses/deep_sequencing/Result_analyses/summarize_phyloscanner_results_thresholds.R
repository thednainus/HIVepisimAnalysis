#summarize results from phyloscanner
#based on the whether two pairs are considered linked (independent of directness)
library(ape)
library(stringr)
library(lubridate)
library(dplyr)
library(HIVepisimAnalysis)
library(tidyverse)
library(caret)

get_prc_data <- function(df_data){

  thresholds <- seq(from = 0, to = 1, by = 0.01)
  rates <- lapply(thresholds, get_prc_deepseq, df_data)
  rates <- do.call(rbind, rates)

  return(rates)

}

get_prc_deepseq <- function(threshold, df_true, byReplicate = FALSE){

  #browser()

  param <- unique(df_true$param)
  mig <- unique(df_true$mig)
  #perc = unique(df_true$perc)
  #code = unique(df_true$code)
  #sampler = unique(df_true$sampler)
  mig = unique(df_true$mig)

  if(byReplicate == FALSE){
    rep = "merged"
  } else {
    rep = unique(df_true$rep)
  }


  df_true["pred"] <- ifelse(df_true$prop_linked >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  df_true$labels <- as.character(df_true$labels)
  df_true$labels <- as.factor(df_true$labels)

  positives <- sum(df_true$labels == 1)
  negatives <- length(df_true$labels) - positives

  #check whether the data can be used with the function caret:confusionMatrix

  #check_data_cm <- check_data_cm(data = df_true$pred,
  #                               reference = df_true$labels,
  #                               positive = "1",
  #                               mode = "prec_recall")
  check_data_cm <- "ok"
  if(check_data_cm == "ok"){

    #create confusion matrix
    cm <- confusionMatrix(df_true$pred, df_true$labels,
                          positive = "1",
                          mode = "prec_recall")
    #true positive rate
    precision <- cm$byClass[[5]]
    #false positive rate
    recall <- cm$byClass[[6]]

    #browser()


    rates <- data.frame(threshold = threshold, precision = precision,
                        recall = recall,
                        positives = positives,
                        negatives = negatives,
                        param = param, rep = rep,
                        mig = mig)

  }

  if(check_data_cm == "not_ok"){


    rates <- data.frame(threshold = threshold, precision = precision,
                        recall = recall,
                        positives = positives,
                        negatives = negatives,
                        param = param, rep = rep,
                        mig = mig)

  }

  return(rates)
}



#filtering by W > 0.01
all_phyloscanner_results_W0.01 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_phyloscanner_results_test_W0.01_final_TN_TP_sampleObserved.RDS")
all_phyloscanner_results_W0.01["group"] <- paste(all_phyloscanner_results_W0.01$param,
                                                 all_phyloscanner_results_W0.01$mig,
                                                 sep = "_")
all_phyloscanner_results_W0.01["reps_groups"] <- paste(all_phyloscanner_results_W0.01$mig,
                                                       all_phyloscanner_results_W0.01$param,
                                                       all_phyloscanner_results_W0.01$rep,
                                                       sep = "_")

all_phyloscanner_results_W0.01["labels"] <- ifelse(all_phyloscanner_results_W0.01$real_pair == "yes", 1, 0)
all_phyloscanner_results_W0.01["tag"] <- paste(all_phyloscanner_results_W0.01$param,
                                               all_phyloscanner_results_W0.01$mig, sep = "_")


prc_teste <- all_phyloscanner_results_W0.01 %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })

prc_all_w0.01 <- do.call(rbind, prc_teste)
saveRDS(prc_all_w0.01, "prc_all_deepseq_w0.01.RDS")




