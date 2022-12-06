# summary of TP (true positive), FP (false positive), TN (true negative) and
# FN (false negative) for consensus sequences.
# Here we assumed that a true transmission pair (independent of who infected whom),
# was a pair that showed infector probability >= 90% or infector probability >= 90%

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


# Sampler 1 250mig threshold 90%----
#read data
#migration = 250
TrueTrees_250mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_250mig_threshold0.90.RDS")
ML1000bp_250mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_1000bp_250mig_threshold0.90.RDS")
ML10000bp_250mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_10000bp_250mig_threshold0.90.RDS")


s1_250mig_90 <- rbind(TrueTrees_250mig_s1_90, ML1000bp_250mig_s1_90, ML10000bp_250mig_s1_90)

s1_250mig_90_obs <- add_abserved_values(s1_250mig_90)
s1_250mig_90_obs["groupby"] <- paste(s1_250mig_90_obs$param_perc,
                                     s1_250mig_90_obs$code,
                                     sep = "_")

#total of who infected whom
wiw_s1_90 <- s1_250mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s1_ordered_90 <- wiw_s1_90[order(wiw_s1_90$perc,
                                     wiw_s1_90$param,
                                     wiw_s1_90$code),]


results_consensus_s1_90_250mig <- s1_250mig_90_obs %>%
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

results_consensus_s1_90_250mig["total"] <- rowSums(results_consensus_s1_90_250mig[,8:11])

results_consensus_s1_90_250mig["TP_perc"] <- results_consensus_s1_90_250mig$total_TP/results_consensus_s1_90_250mig$total
results_consensus_s1_90_250mig["FP_perc"] <- results_consensus_s1_90_250mig$total_FP/results_consensus_s1_90_250mig$total
results_consensus_s1_90_250mig["TN_perc"] <- results_consensus_s1_90_250mig$total_TN/results_consensus_s1_90_250mig$total
results_consensus_s1_90_250mig["FN_perc"] <- results_consensus_s1_90_250mig$total_FN/results_consensus_s1_90_250mig$total

all_250mig_s1_90data <- results_consensus_s1_90_250mig[order(results_consensus_s1_90_250mig$perc,
                                                             results_consensus_s1_90_250mig$param,
                                                             results_consensus_s1_90_250mig$code),]

all_250mig_s1_90_data2 <- all_250mig_s1_90data[,c(1:12)]
all_250mig_s1_90_data2["total_wiw_within_TP"] <- wiw_s1_ordered_90$total_wiw


# Sampler 2 250mig threshold 90%----
#read data
#migration = 250
TrueTrees_250mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_250mig_threshold0.90.RDS")
ML1000bp_250mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_1000bp_250mig_threshold0.90.RDS")
ML10000bp_250mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_10000bp_250mig_threshold0.90.RDS")


s2_250mig_90 <- rbind(TrueTrees_250mig_s2_90, ML1000bp_250mig_s2_90, ML10000bp_250mig_s2_90)

s2_250mig_90_obs <- add_abserved_values(s2_250mig_90)
s2_250mig_90_obs["groupby"] <- paste(s2_250mig_90_obs$param_perc,
                                     s2_250mig_90_obs$code,
                                     sep = "_")

#total of who infected whom
wiw_s2_90 <- s2_250mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s2_ordered_90 <- wiw_s2_90[order(wiw_s2_90$perc,
                                     wiw_s2_90$param,
                                     wiw_s2_90$code),]


results_consensus_s2_90_250mig <- s2_250mig_90_obs %>%
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

results_consensus_s2_90_250mig["total"] <- rowSums(results_consensus_s2_90_250mig[,8:11])

results_consensus_s2_90_250mig["TP_perc"] <- results_consensus_s2_90_250mig$total_TP/results_consensus_s2_90_250mig$total
results_consensus_s2_90_250mig["FP_perc"] <- results_consensus_s2_90_250mig$total_FP/results_consensus_s2_90_250mig$total
results_consensus_s2_90_250mig["TN_perc"] <- results_consensus_s2_90_250mig$total_TN/results_consensus_s2_90_250mig$total
results_consensus_s2_90_250mig["FN_perc"] <- results_consensus_s2_90_250mig$total_FN/results_consensus_s2_90_250mig$total

all_250mig_s2_90data <- results_consensus_s2_90_250mig[order(results_consensus_s2_90_250mig$perc,
                                                             results_consensus_s2_90_250mig$param,
                                                             results_consensus_s2_90_250mig$code),]

all_250mig_s2_90_data2 <- all_250mig_s2_90data[,c(1:12)]
all_250mig_s2_90_data2["total_wiw_within_TP"] <- wiw_s2_ordered_90$total_wiw












# Sampler 1 500mig threshold 90%----
#read data
#migration = 500
TrueTrees_500mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_500mig_threshold0.90.RDS")
ML1000bp_500mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_1000bp_500mig_threshold0.90.RDS")
ML10000bp_500mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_10000bp_500mig_threshold0.90.RDS")


s1_500mig_90 <- rbind(TrueTrees_500mig_s1_90, ML1000bp_500mig_s1_90, ML10000bp_500mig_s1_90)

s1_500mig_90_obs <- add_abserved_values(s1_500mig_90)
s1_500mig_90_obs["groupby"] <- paste(s1_500mig_90_obs$param_perc,
                                      s1_500mig_90_obs$code,
                                      sep = "_")

#total of who infected whom
wiw_s1_90 <- s1_500mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s1_ordered_90 <- wiw_s1_90[order(wiw_s1_90$perc,
                                     wiw_s1_90$param,
                                     wiw_s1_90$code),]


results_consensus_s1_90_500mig <- s1_500mig_90_obs %>%
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

results_consensus_s1_90_500mig["total"] <- rowSums(results_consensus_s1_90_500mig[,8:11])

results_consensus_s1_90_500mig["TP_perc"] <- results_consensus_s1_90_500mig$total_TP/results_consensus_s1_90_500mig$total
results_consensus_s1_90_500mig["FP_perc"] <- results_consensus_s1_90_500mig$total_FP/results_consensus_s1_90_500mig$total
results_consensus_s1_90_500mig["TN_perc"] <- results_consensus_s1_90_500mig$total_TN/results_consensus_s1_90_500mig$total
results_consensus_s1_90_500mig["FN_perc"] <- results_consensus_s1_90_500mig$total_FN/results_consensus_s1_90_500mig$total

all_500mig_s1_90data <- results_consensus_s1_90_500mig[order(results_consensus_s1_90_500mig$perc,
                                                             results_consensus_s1_90_500mig$param,
                                                             results_consensus_s1_90_500mig$code),]

all_500mig_s1_90_data2 <- all_500mig_s1_90data[,c(1:12)]
all_500mig_s1_90_data2["total_wiw_within_TP"] <- wiw_s1_ordered_90$total_wiw


# Sampler 2 500mig threshold 90%----
#read data
#migration = 500
TrueTrees_500mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_500mig_threshold0.90.RDS")
ML1000bp_500mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_1000bp_500mig_threshold0.90.RDS")
ML10000bp_500mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_10000bp_500mig_threshold0.90.RDS")


s2_500mig_90 <- rbind(TrueTrees_500mig_s2_90, ML1000bp_500mig_s2_90, ML10000bp_500mig_s2_90)

s2_500mig_90_obs <- add_abserved_values(s2_500mig_90)
s2_500mig_90_obs["groupby"] <- paste(s2_500mig_90_obs$param_perc,
                                     s2_500mig_90_obs$code,
                                     sep = "_")

#total of who infected whom
wiw_s2_90 <- s2_500mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s2_ordered_90 <- wiw_s2_90[order(wiw_s2_90$perc,
                                     wiw_s2_90$param,
                                     wiw_s2_90$code),]


results_consensus_s2_90_500mig <- s2_500mig_90_obs %>%
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

results_consensus_s2_90_500mig["total"] <- rowSums(results_consensus_s2_90_500mig[,8:11])

results_consensus_s2_90_500mig["TP_perc"] <- results_consensus_s2_90_500mig$total_TP/results_consensus_s2_90_500mig$total
results_consensus_s2_90_500mig["FP_perc"] <- results_consensus_s2_90_500mig$total_FP/results_consensus_s2_90_500mig$total
results_consensus_s2_90_500mig["TN_perc"] <- results_consensus_s2_90_500mig$total_TN/results_consensus_s2_90_500mig$total
results_consensus_s2_90_500mig["FN_perc"] <- results_consensus_s2_90_500mig$total_FN/results_consensus_s2_90_500mig$total

all_500mig_s2_90data <- results_consensus_s2_90_500mig[order(results_consensus_s2_90_500mig$perc,
                                                             results_consensus_s2_90_500mig$param,
                                                             results_consensus_s2_90_500mig$code),]

all_500mig_s2_90_data2 <- all_500mig_s2_90data[,c(1:12)]
all_500mig_s2_90_data2["total_wiw_within_TP"] <- wiw_s2_ordered_90$total_wiw





# Sampler 1 750mig threshold 90%----
#read data
#migration = 750
TrueTrees_750mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_750mig_threshold0.90.RDS")
ML1000bp_750mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_1000bp_750mig_threshold0.90.RDS")
ML10000bp_750mig_s1_90 <- readRDS("Analyses/consensus_sequences/Results/sampler1/threshold_0.90/all_data_s1_10000bp_750mig_threshold0.90.RDS")


s1_750mig_90 <- rbind(TrueTrees_750mig_s1_90, ML1000bp_750mig_s1_90, ML10000bp_750mig_s1_90)

s1_750mig_90_obs <- add_abserved_values(s1_750mig_90)
s1_750mig_90_obs["groupby"] <- paste(s1_750mig_90_obs$param_perc,
                                     s1_750mig_90_obs$code,
                                     sep = "_")

#total of who infected whom
wiw_s1_90 <- s1_750mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s1_ordered_90 <- wiw_s1_90[order(wiw_s1_90$perc,
                                     wiw_s1_90$param,
                                     wiw_s1_90$code),]


results_consensus_s1_90_750mig <- s1_750mig_90_obs %>%
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

results_consensus_s1_90_750mig["total"] <- rowSums(results_consensus_s1_90_750mig[,8:11])

results_consensus_s1_90_750mig["TP_perc"] <- results_consensus_s1_90_750mig$total_TP/results_consensus_s1_90_750mig$total
results_consensus_s1_90_750mig["FP_perc"] <- results_consensus_s1_90_750mig$total_FP/results_consensus_s1_90_750mig$total
results_consensus_s1_90_750mig["TN_perc"] <- results_consensus_s1_90_750mig$total_TN/results_consensus_s1_90_750mig$total
results_consensus_s1_90_750mig["FN_perc"] <- results_consensus_s1_90_750mig$total_FN/results_consensus_s1_90_750mig$total

all_750mig_s1_90data <- results_consensus_s1_90_750mig[order(results_consensus_s1_90_750mig$perc,
                                                             results_consensus_s1_90_750mig$param,
                                                             results_consensus_s1_90_750mig$code),]

all_750mig_s1_90_data2 <- all_750mig_s1_90data[,c(1:12)]
all_750mig_s1_90_data2["total_wiw_within_TP"] <- wiw_s1_ordered_90$total_wiw


# Sampler 2 750mig threshold 90%----
#read data
#migration = 750
TrueTrees_750mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_750mig_threshold0.90.RDS")
ML1000bp_750mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_1000bp_750mig_threshold0.90.RDS")
ML10000bp_750mig_s2_90 <- readRDS("Analyses/consensus_sequences/Results/sampler2/threshold_0.90/all_data_s2_10000bp_750mig_threshold0.90.RDS")


s2_750mig_90 <- rbind(TrueTrees_750mig_s2_90, ML1000bp_750mig_s2_90, ML10000bp_750mig_s2_90)

s2_750mig_90_obs <- add_abserved_values(s2_750mig_90)
s2_750mig_90_obs["groupby"] <- paste(s2_750mig_90_obs$param_perc,
                                     s2_750mig_90_obs$code,
                                     sep = "_")

#total of who infected whom
wiw_s2_90 <- s2_750mig_90_obs %>%
  group_by(groupby) %>%
  mutate(total_wiw = sum(labels == 1 & linked == "yes" & real_pair == "yes")) %>%
  select(groupby, param, perc, code, total_wiw) %>%
  distinct()

wiw_s2_ordered_90 <- wiw_s2_90[order(wiw_s2_90$perc,
                                     wiw_s2_90$param,
                                     wiw_s2_90$code),]


results_consensus_s2_90_750mig <- s2_750mig_90_obs %>%
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

results_consensus_s2_90_750mig["total"] <- rowSums(results_consensus_s2_90_750mig[,8:11])

results_consensus_s2_90_750mig["TP_perc"] <- results_consensus_s2_90_750mig$total_TP/results_consensus_s2_90_750mig$total
results_consensus_s2_90_750mig["FP_perc"] <- results_consensus_s2_90_750mig$total_FP/results_consensus_s2_90_750mig$total
results_consensus_s2_90_750mig["TN_perc"] <- results_consensus_s2_90_750mig$total_TN/results_consensus_s2_90_750mig$total
results_consensus_s2_90_750mig["FN_perc"] <- results_consensus_s2_90_750mig$total_FN/results_consensus_s2_90_750mig$total

all_750mig_s2_90data <- results_consensus_s2_90_750mig[order(results_consensus_s2_90_750mig$perc,
                                                             results_consensus_s2_90_750mig$param,
                                                             results_consensus_s2_90_750mig$code),]

all_750mig_s2_90_data2 <- all_750mig_s2_90data[,c(1:12)]
all_750mig_s2_90_data2["total_wiw_within_TP"] <- wiw_s2_ordered_90$total_wiw

