#get total tips analysed in each population dept (5%, 10%, to 90%)

library(dplyr)

#sampler 1: 1,000bp
mig250_s1_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s1_1000bp_250mig_tips.RDS")
mig250_s1_1000bp["mig"] <- "1/0.68"
mig250_s1_1000bp["perc"] <- as.factor(as.numeric(mig250_s1_1000bp$perc) * 100)
mig250_s1_1000bp["param"] <- ifelse(mig250_s1_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")
mig500_s1_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s1_1000bp_500mig_tips.RDS")
mig500_s1_1000bp["mig"] <- "1/1.37"
mig500_s1_1000bp["perc"] <- as.factor(as.numeric(mig500_s1_1000bp$perc) * 100)
mig500_s1_1000bp["param"] <- ifelse(mig500_s1_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")
mig750_s1_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s1_1000bp_750mig_tips.RDS")
mig750_s1_1000bp["mig"] <- "1/2.05"
mig750_s1_1000bp["perc"] <- as.factor(as.numeric(mig750_s1_1000bp$perc) * 100)
mig750_s1_1000bp["param"] <- ifelse(mig750_s1_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")
#all sampler 1
#all_s1_1000bp <- rbind(mig250_s1_1000bp, mig500_s1_1000bp, mig750_s1_1000bp)

#sampler 2: 1,000bp
mig250_s2_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s2_1000bp_250mig_tips.RDS")
mig250_s2_1000bp["mig"] <- "1/0.68"
mig250_s2_1000bp["perc"] <- as.factor(as.numeric(mig250_s2_1000bp$perc) * 100)
mig250_s2_1000bp["param"] <- ifelse(mig250_s2_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")
mig500_s2_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s2_1000bp_500mig_tips.RDS")
mig500_s2_1000bp["mig"] <- "1/1.37"
mig500_s2_1000bp["perc"] <- as.factor(as.numeric(mig500_s2_1000bp$perc) * 100)
mig500_s2_1000bp["param"] <- ifelse(mig500_s2_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")
mig750_s2_1000bp <- readRDS("~/OneDrive - Imperial College London/HIV_SanDiego/tip_numbers/all_data_s2_1000bp_750mig_tips.RDS")
mig750_s2_1000bp["mig"] <- "1/2.05"
mig750_s2_1000bp["perc"] <- as.factor(as.numeric(mig750_s2_1000bp$perc) * 100)
mig750_s2_1000bp["param"] <- ifelse(mig750_s2_1000bp$param == "1067",
                                    "Combination 1", "Combination 2")

#all sampler 2
#all_s2_1000bp <- rbind(mig250_s2_1000bp, mig500_s2_1000bp, mig750_s2_1000bp)




# sampler 1
summary_tips_s1_250mig <- mig250_s1_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                           sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()

summary_tips_s1_500mig <- mig500_s1_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                        sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()

summary_tips_s1_750mig <- mig750_s1_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                        sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()


# sampler 2
summary_tips_s2_250mig <- mig250_s2_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                        sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()

summary_tips_s2_500mig <- mig500_s2_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                        sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()

summary_tips_s2_750mig <- mig750_s2_1000bp %>% group_by(param_perc) %>%
  mutate(mean = summary(tree_size_region)[[4]],
         median = summary(tree_size_region)[[3]],
         lower_bound = summary(tree_size_region)[[2]],
         upper_bound = summary(tree_size_region)[[5]]) %>%
  select(param_perc, param, perc, mig, median, lower_bound, upper_bound) %>%
  mutate(median = paste(median, " (CI = ", lower_bound, "--", upper_bound, ")",
                        sep = ""))  %>%
  select(param_perc, param, perc, mig, median) %>%
  distinct()





library(xtable)

print(xtable(summary_tips_s1_250mig[,2:5]), include.rownames=FALSE)
print(xtable(summary_tips_s1_500mig[,2:5]), include.rownames=FALSE)
print(xtable(summary_tips_s1_750mig[,2:5]), include.rownames=FALSE)

print(xtable(summary_tips_s2_250mig[,2:5]), include.rownames=FALSE)
print(xtable(summary_tips_s2_500mig[,2:5]), include.rownames=FALSE)
print(xtable(summary_tips_s2_750mig[,2:5]), include.rownames=FALSE)
