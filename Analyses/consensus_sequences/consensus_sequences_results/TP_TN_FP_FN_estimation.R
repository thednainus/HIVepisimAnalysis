# summary of TP (true positive), FP (false positive), TN (true negative) and
# FN (false negative) for consensus sequences.
# Here we assumed that a true transmission pair (independent of who infected whom),
# was a pair that showed infector probability >= 80%

library(HIVepisimAnalysis)
library(dplyr)


#function to add observed values of TP, FP, TN, FN

add_abserved_values <- function(df1){

  df1["observed"] <- 0

  df1$observed[df1$linked == "yes" & df1$real_pair == "yes"] <- "TP"
  df1$observed[df1$linked == "no" & df1$real_pair == "yes"] <- "FN"
  df1$observed[df1$linked == "yes" & df1$real_pair == "no"] <- "FP"
  df1$observed[df1$linked == "no" & df1$real_pair == "no"] <- "TN"


  return(df1)

}

#get the total number of who infected whom within the True positives
get_total_wiw <- function(df1){

  toatl_wiw <- sum(df1$labels == 1 & df1$linked == "yes")
  if(df1$labels == 1 & df1$linked == "yes" & df1$real_pair == "yes"){




  }


}


# Sampler 1 ----
#read data
TrueTrees_500mig_1 <- readRDS("consensus_sequences_results/all_data_s1_500mig.RDS")
TrueTrees_500mig_2 <- readRDS("consensus_sequences_results/all_data_s1_500mig_perc_0.05.RDS")
ML1000bp_500mig_1 <- readRDS("consensus_sequences_results/all_data_s1_1000bp_500mig.RDS")
ML1000bp_500mig_2 <- readRDS("consensus_sequences_results/all_data_s1_1000bp_500mig_perc_0.05.RDS")
ML10000bp_500mig_1 <- readRDS("consensus_sequences_results/all_data_s1_10000bp_500mig.RDS")
ML10000bp_500mig_2 <- readRDS("consensus_sequences_results/all_data_s1_10000bp_500mig_perc_0.05.RDS")

s1_500mig <- rbind(TrueTrees_500mig_1, TrueTrees_500mig_2,
                   ML1000bp_500mig_1, ML1000bp_500mig_2,
                   ML10000bp_500mig_1, ML10000bp_500mig_2)

s1_500mig_obs <- add_abserved_values(s1_500mig)
s1_500mig_obs["groupby"] <- paste(s1_500mig_obs$param_perc,
                                  s1_500mig_obs$code,
                                  sep = "_")

#total of who infected whom
wiw <- s1_500mig_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_ordered <- wiw[order(wiw$perc,
                         wiw$param,
                         wiw$code),]


results_consensus_s1_500mig <- s1_500mig_obs %>%
  group_by(groupby) %>%
  mutate(precision = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FP")),
         sensitivity = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FN")),
         specificity = sum(observed == "TN")/(sum(observed == "TN") + sum(observed == "FP")),
         total_TP = sum(observed == "TP"),
         total_FP = sum(observed == "FP"),
         total_TN = sum(observed == "TN"),
         total_FN = sum(observed == "FN")) %>%
  select(groupby, param, perc, code, precision, sensitivity, specificity,
         total_TP, total_FP, total_TN, total_FN) %>%
  distinct()

results_consensus_s1_500mig["total"] <- rowSums(results_consensus_s1_500mig[,8:11])

results_consensus_s1_500mig["TP_perc"] <- results_consensus_s1_500mig$total_TP/results_consensus_s1_500mig$total
results_consensus_s1_500mig["FP_perc"] <- results_consensus_s1_500mig$total_FP/results_consensus_s1_500mig$total
results_consensus_s1_500mig["TN_perc"] <- results_consensus_s1_500mig$total_TN/results_consensus_s1_500mig$total
results_consensus_s1_500mig["FN_perc"] <- results_consensus_s1_500mig$total_FN/results_consensus_s1_500mig$total

all_500mig_data <- results_consensus_s1_500mig[order(results_consensus_s1_500mig$perc,
                                                     results_consensus_s1_500mig$param,
                                                     results_consensus_s1_500mig$code),]

all_500mig_data2 <- all_500mig_data[,c(1:12)]
all_500mig_data2["total_wiw_within_TP"] <- wiw_ordered$total_wiw


# Sampler 2 ----
#read data
TrueTrees_s2_500mig_1 <- readRDS("consensus_sequences_results/all_data_s2_500mig.RDS")
TrueTrees_s2_500mig_2 <- readRDS("consensus_sequences_results/all_data_s2_500mig_perc_0.05.RDS")
ML1000bp_s2_500mig_1 <- readRDS("consensus_sequences_results/all_data_s2_1000bp_500mig.RDS")
ML1000bp_s2_500mig_2 <- readRDS("consensus_sequences_results/all_data_s2_1000bp_500mig_perc_0.05.RDS")
ML10000bp_s2_500mig_1 <- readRDS("consensus_sequences_results/all_data_s2_10000bp_500mig.RDS")
ML10000bp_s2_500mig_2 <- readRDS("consensus_sequences_results/all_data_s2_10000bp_500mig_perc_0.05.RDS")

s2_500mig <- rbind(TrueTrees_s2_500mig_1, TrueTrees_s2_500mig_2,
                   ML1000bp_s2_500mig_1, ML1000bp_s2_500mig_2,
                   ML10000bp_s2_500mig_1, ML10000bp_s2_500mig_2)

s2_500mig_obs <- add_abserved_values(s2_500mig)
s2_500mig_obs["groupby"] <- paste(s2_500mig_obs$param_perc,
                                  s2_500mig_obs$code,
                                  sep = "_")

#total of who infected whom
wiw_s2 <- s2_500mig_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_ordered_s2 <- wiw_s2[order(wiw_s2$perc,
                               wiw_s2$param,
                               wiw_s2$code),]




results_consensus_s2_500mig <- s2_500mig_obs %>%
  group_by(groupby) %>%
  mutate(precision = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FP")),
         sensitivity = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FN")),
         specificity = sum(observed == "TN")/(sum(observed == "TN") + sum(observed == "FP")),
         total_TP = sum(observed == "TP"),
         total_FP = sum(observed == "FP"),
         total_TN = sum(observed == "TN"),
         total_FN = sum(observed == "FN")) %>%
  select(groupby, param, perc, code, precision, sensitivity, specificity,
         total_TP, total_FP, total_TN, total_FN) %>%
  distinct()

results_consensus_s2_500mig["total"] <- rowSums(results_consensus_s2_500mig[,8:11])

results_consensus_s2_500mig["TP_perc"] <- results_consensus_s2_500mig$total_TP/results_consensus_s2_500mig$total
results_consensus_s2_500mig["FP_perc"] <- results_consensus_s2_500mig$total_FP/results_consensus_s2_500mig$total
results_consensus_s2_500mig["TN_perc"] <- results_consensus_s2_500mig$total_TN/results_consensus_s2_500mig$total
results_consensus_s2_500mig["FN_perc"] <- results_consensus_s2_500mig$total_FN/results_consensus_s2_500mig$total

all_500mig_s2_data <- results_consensus_s2_500mig[order(results_consensus_s2_500mig$perc,
                                                        results_consensus_s2_500mig$param,
                                                        results_consensus_s2_500mig$code),]


all_500mig_s2_data2 <- all_500mig_s2_data[,c(1:12)]
all_500mig_s2_data2["total_wiw_within_TP"] <- wiw_ordered_s2$total_wiw

all_500mig_s2_data2$total_wiw_within_TP/all_500mig_s2_data2$total_TP
