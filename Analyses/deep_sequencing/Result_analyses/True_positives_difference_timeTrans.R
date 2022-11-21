#checking the true positives

library(dplyr)

#read pairs in which filtering was W >= 80%

W80 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_phyloscanner_results_test_W80_final_TN_TP_sampleObserved.RDS")
W80["param_mig"] <- paste(W80$param, W80$mig, sep = "_")
W80["byGroup"] <- c(1:nrow(W80))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)

TP_W80 <- subset(W80, observed == "TP")
TP_W80["difference_donor_tt"] <- (TP_W80$st_donor_ID - TP_W80$time_trans)*12
TP_W80["difference_recip_tt"] <- (TP_W80$st_recip_ID - TP_W80$time_trans)*12
TP_W80["difference_donor_recip"] <- (TP_W80$st_donor_ID - TP_W80$st_recip_ID)*12

TP_W80$difference_donor_recip <- abs(TP_W80$difference_donor_recip)




TP_W80 %>%
  group_by(param_mig) %>%
  mutate(mean = mean(difference_donor_recip),
         median = median(difference_donor_recip),
         Min = summary(difference_donor_recip)[1],
         Max = summary(difference_donor_recip)[6],
         quantile1 = summary(difference_donor_recip)[2],
         median1 = summary(difference_donor_recip)[3],
         quantile3 = summary(difference_donor_recip)[5],
         total_rows = n()) %>%
  select(param_mig, mean, median, median1, total_rows, Min, Max, quantile1, quantile3) %>%
  distinct()

