#get the percentage of true positives that is a transmission pair
#and that phyloscanner correctly identified who infected whom

library(dplyr)
library(HIVepisimAnalysis)

#filtering by W > 0.80
all_phyloscanner_results_W80 <- readRDS("all_phyloscanner_results_test_W80_final.RDS")
all_phyloscanner_results_W80["group"] <- paste(all_phyloscanner_results_W80$param,
                                               all_phyloscanner_results_W80$mig,
                                               sep = "_")
all_phyloscanner_results_W80["reps_groups"] <- paste(all_phyloscanner_results_W80$mig,
                                                     all_phyloscanner_results_W80$param,
                                                     all_phyloscanner_results_W80$rep,
                                                     sep = "_")


all_phyloscanner_results_W80 %>%
  group_by(group) %>%
  mutate(total_correct = sum(observed == "TP" &
                          real_pair == "yes" &
                          directness == "no")) %>%
  select(group, total_correct) %>%
  distinct()



TP_correct_W80 <- all_phyloscanner_results_W80 %>%
  group_by(group) %>%
  group_modify(~ wiw(.x))

TP_correct_W80["percentage"] <- TP_correct_W80$total_correct/TP_correct_W80$total_TP



#filtering by W > 0.01
all_phyloscanner_results_W0.01 <- readRDS("all_phyloscanner_results_test_W0.01_final.RDS")
all_phyloscanner_results_W0.01["group"] <- paste(all_phyloscanner_results_W0.01$param,
                                                 all_phyloscanner_results_W0.01$mig,
                                                 sep = "_")
all_phyloscanner_results_W0.01["reps_groups"] <- paste(all_phyloscanner_results_W0.01$mig,
                                                       all_phyloscanner_results_W0.01$param,
                                                       all_phyloscanner_results_W0.01$rep,
                                                       sep = "_")



TP_correct_W0.01 <- all_phyloscanner_results_W0.01 %>%
  group_by(group) %>%
  group_modify(~ wiw(.x))

TP_correct_W0.01["percentage"] <- TP_correct_W0.01$total_correct/TP_correct_W0.01$total_TP
