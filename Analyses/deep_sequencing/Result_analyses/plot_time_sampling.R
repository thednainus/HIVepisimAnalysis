#plot time of sampling
library(ggplot2)

#read data
phyloscanner_summary_results <- readRDS("all_phyloscanner_results_test1.RDS")

#subset_data
#difference is in year. I multiply by 12 to get the units in month for the differences
# between sample time of donor or recipient and time of transmission
results <- subset(phyloscanner_summary_results, trans == "true")
results["difference_donor_tt"] <- (results$st_donor_ID - results$time_trans)*12
results["difference_recip_tt"] <- (results$st_recip_ID - results$time_trans)*12
results["difference_donor_recip"] <- (results$st_donor_ID - results$st_recip_ID)*12

results["pair_number"] <- c(1:nrow(results))
results$pair_number <- as.character(results$pair_number)
results$pair_number <- as.factor(results$pair_number)


results1 <- subset(results, mig == "250migrants" & param == "2348")

quartz()
ggplot(results1, aes(x = difference_recip_tt, y = pair_number)) +
  geom_point(aes(col = direction)) +
  ggtitle("250 mig: param 1067: difference_recip_tt (months) ")

quartz()
ggplot(results1, aes(x = difference_donor_tt, y = pair_number)) +
  geom_point(aes(col = direction)) +
  ggtitle("250 mig: param 1067: difference_donor_tt (months) ")

quartz()
ggplot(results1, aes(x = difference_donor_recip, y = pair_number)) +
  geom_point(aes(col = direction)) +
  ggtitle("250 mig: param 1067: difference_donor_recip (months) ")




