# Summarize results for consensus sequences
# It will read the labels results (1 and 0 to determine if correctly detected
# a transmission pair or not)

library(tidyverse)
library(HIVepisimAnalysis)
library(xtable)

# Sampler 1----
#for true trees
mig250_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_250mig.RDS")
mig250_true["mig"] <- "1/0.68"
mig250_true["perc"] <- as.factor(as.numeric(mig250_true$perc) * 100)
mig250_true["param"] <- ifelse(mig250_true$param == "1067",
                               "Combination 1", "Combination 2")
mig250_true["tag"] <- paste(mig250_true$param_perc, mig250_true$mig, sep = "_")

mig500_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_500mig.RDS")
mig500_true["mig"] <- "1/1.37"
mig500_true["perc"] <- as.factor(as.numeric(mig500_true$perc) * 100)
mig500_true["param"] <- ifelse(mig500_true$param == "1067",
                               "Combination 1", "Combination 2")
mig500_true["tag"] <- paste(mig500_true$param_perc, mig500_true$mig, sep = "_")

mig750_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_750mig.RDS")
mig750_true["mig"] <- "1/2.05"
mig750_true["perc"] <- as.factor(as.numeric(mig750_true$perc) * 100)
mig750_true["param"] <- ifelse(mig750_true$param == "1067",
                               "Combination 1", "Combination 2")
mig750_true["tag"] <- paste(mig750_true$param_perc, mig750_true$mig, sep = "_")



#ML 1,000bp
mig250_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_250mig.RDS")
mig250_1000bp["mig"] <- "1/0.68"
mig250_1000bp["perc"] <- as.factor(as.numeric(mig250_1000bp$perc) * 100)
mig250_1000bp["param"] <- ifelse(mig250_1000bp$param == "1067",
                               "Combination 1", "Combination 2")
mig250_1000bp["tag"] <- paste(mig250_1000bp$param_perc, mig250_1000bp$mig, sep = "_")

mig500_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_500mig.RDS")
mig500_1000bp["mig"] <- "1/1.37"
mig500_1000bp["perc"] <- as.factor(as.numeric(mig500_1000bp$perc) * 100)
mig500_1000bp["param"] <- ifelse(mig500_1000bp$param == "1067",
                               "Combination 1", "Combination 2")
mig500_1000bp["tag"] <- paste(mig500_1000bp$param_perc, mig500_1000bp$mig, sep = "_")

mig750_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_750mig.RDS")
mig750_1000bp["mig"] <- "1/2.05"
mig750_1000bp["perc"] <- as.factor(as.numeric(mig750_1000bp$perc) * 100)
mig750_1000bp["param"] <- ifelse(mig750_1000bp$param == "1067",
                               "Combination 1", "Combination 2")
mig750_1000bp["tag"] <- paste(mig750_1000bp$param_perc, mig750_1000bp$mig, sep = "_")


#ML 10,000bp
mig250_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_250mig.RDS")
mig250_10000bp["mig"] <- "1/0.68"
mig250_10000bp["perc"] <- as.factor(as.numeric(mig250_10000bp$perc) * 100)
mig250_10000bp["param"] <- ifelse(mig250_10000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig250_10000bp["tag"] <- paste(mig250_10000bp$param_perc, mig250_10000bp$mig, sep = "_")


mig500_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_500mig.RDS")
mig500_10000bp["mig"] <- "1/1.37"
mig500_10000bp["perc"] <- as.factor(as.numeric(mig500_10000bp$perc) * 100)
mig500_10000bp["param"] <- ifelse(mig500_10000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig500_10000bp["tag"] <- paste(mig500_10000bp$param_perc, mig500_10000bp$mig, sep = "_")


mig750_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_750mig.RDS")
mig750_10000bp["mig"] <- "1/2.05"
mig750_10000bp["perc"] <- as.factor(as.numeric(mig750_10000bp$perc) * 100)
mig750_10000bp["param"] <- ifelse(mig750_10000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig750_10000bp["tag"] <- paste(mig750_10000bp$param_perc, mig750_10000bp$mig, sep = "_")



mig250_true_s1_totalPairs <- mig250_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250_1000bp_s1_totalPairs <- mig250_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250_10000bp_s1_totalPairs <- mig250_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250 <- cbind(mig250_true_s1_totalPairs,
                mig250_1000bp_s1_totalPairs$total,
                mig250_10000bp_s1_totalPairs$total)

names(mig250)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig250[,c(3:5,8:10)]), include.rownames=FALSE)





mig500_true_s1_totalPairs <- mig500_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig500_1000bp_s1_totalPairs <- mig500_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig500_10000bp_s1_totalPairs <- mig500_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()





mig500 <- cbind(mig500_true_s1_totalPairs,
                mig500_1000bp_s1_totalPairs$total,
                mig500_10000bp_s1_totalPairs$total)

names(mig500)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig500[,c(3:5,8:10)]), include.rownames=FALSE)




mig750_true_s1_totalPairs <- mig750_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig750_1000bp_s1_totalPairs <- mig750_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig750_10000bp_s1_totalPairs <- mig750_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()





mig750 <- cbind(mig750_true_s1_totalPairs,
                mig750_1000bp_s1_totalPairs$total,
                mig750_10000bp_s1_totalPairs$total)

names(mig750)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig750[,c(3:5,8:10)]), include.rownames=FALSE)









# sampler2 ----
#for true trees
mig250_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_250mig.RDS")
mig250_true["mig"] <- "1/0.68"
mig250_true["perc"] <- as.factor(as.numeric(mig250_true$perc) * 100)
mig250_true["param"] <- ifelse(mig250_true$param == "1067",
                               "Combination 1", "Combination 2")
mig250_true["tag"] <- paste(mig250_true$param_perc, mig250_true$mig, sep = "_")

mig500_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_500mig.RDS")
mig500_true["mig"] <- "1/1.37"
mig500_true["perc"] <- as.factor(as.numeric(mig500_true$perc) * 100)
mig500_true["param"] <- ifelse(mig500_true$param == "1067",
                               "Combination 1", "Combination 2")
mig500_true["tag"] <- paste(mig500_true$param_perc, mig500_true$mig, sep = "_")

mig750_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_750mig.RDS")
mig750_true["mig"] <- "1/2.05"
mig750_true["perc"] <- as.factor(as.numeric(mig750_true$perc) * 100)
mig750_true["param"] <- ifelse(mig750_true$param == "1067",
                               "Combination 1", "Combination 2")
mig750_true["tag"] <- paste(mig750_true$param_perc, mig750_true$mig, sep = "_")



#ML 1,000bp
mig250_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_250mig.RDS")
mig250_1000bp["mig"] <- "1/0.68"
mig250_1000bp["perc"] <- as.factor(as.numeric(mig250_1000bp$perc) * 100)
mig250_1000bp["param"] <- ifelse(mig250_1000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig250_1000bp["tag"] <- paste(mig250_1000bp$param_perc, mig250_1000bp$mig, sep = "_")

mig500_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_500mig.RDS")
mig500_1000bp["mig"] <- "1/1.37"
mig500_1000bp["perc"] <- as.factor(as.numeric(mig500_1000bp$perc) * 100)
mig500_1000bp["param"] <- ifelse(mig500_1000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig500_1000bp["tag"] <- paste(mig500_1000bp$param_perc, mig500_1000bp$mig, sep = "_")

mig750_1000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_750mig.RDS")
mig750_1000bp["mig"] <- "1/2.05"
mig750_1000bp["perc"] <- as.factor(as.numeric(mig750_1000bp$perc) * 100)
mig750_1000bp["param"] <- ifelse(mig750_1000bp$param == "1067",
                                 "Combination 1", "Combination 2")
mig750_1000bp["tag"] <- paste(mig750_1000bp$param_perc, mig750_1000bp$mig, sep = "_")


#ML 10,000bp
mig250_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_250mig.RDS")
mig250_10000bp["mig"] <- "1/0.68"
mig250_10000bp["perc"] <- as.factor(as.numeric(mig250_10000bp$perc) * 100)
mig250_10000bp["param"] <- ifelse(mig250_10000bp$param == "1067",
                                  "Combination 1", "Combination 2")
mig250_10000bp["tag"] <- paste(mig250_10000bp$param_perc, mig250_10000bp$mig, sep = "_")


mig500_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_500mig.RDS")
mig500_10000bp["mig"] <- "1/1.37"
mig500_10000bp["perc"] <- as.factor(as.numeric(mig500_10000bp$perc) * 100)
mig500_10000bp["param"] <- ifelse(mig500_10000bp$param == "1067",
                                  "Combination 1", "Combination 2")
mig500_10000bp["tag"] <- paste(mig500_10000bp$param_perc, mig500_10000bp$mig, sep = "_")


mig750_10000bp <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_750mig.RDS")
mig750_10000bp["mig"] <- "1/2.05"
mig750_10000bp["perc"] <- as.factor(as.numeric(mig750_10000bp$perc) * 100)
mig750_10000bp["param"] <- ifelse(mig750_10000bp$param == "1067",
                                  "Combination 1", "Combination 2")
mig750_10000bp["tag"] <- paste(mig750_10000bp$param_perc, mig750_10000bp$mig, sep = "_")



mig250_true_s2_totalPairs <- mig250_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250_1000bp_s2_totalPairs <- mig250_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250_10000bp_s2_totalPairs <- mig250_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig250 <- cbind(mig250_true_s2_totalPairs,
                mig250_1000bp_s2_totalPairs$total,
                mig250_10000bp_s2_totalPairs$total)

names(mig250)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig250[,c(3:5,8:10)]), include.rownames=FALSE)





mig500_true_s2_totalPairs <- mig500_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig500_1000bp_s2_totalPairs <- mig500_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig500_10000bp_s2_totalPairs <- mig500_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()





mig500 <- cbind(mig500_true_s2_totalPairs,
                mig500_1000bp_s2_totalPairs$total,
                mig500_10000bp_s2_totalPairs$total)

names(mig500)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig500[,c(3:5,8:10)]), include.rownames=FALSE)




mig750_true_s2_totalPairs <- mig750_true %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig750_1000bp_s2_totalPairs <- mig750_1000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()

mig750_10000bp_s2_totalPairs <- mig750_10000bp %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_perc, tag, total) %>%
  distinct()





mig750 <- cbind(mig750_true_s2_totalPairs,
                mig750_1000bp_s2_totalPairs$total,
                mig750_10000bp_s2_totalPairs$total)

names(mig750)[8:10] <- c("True trees", "ML 1,000bp", "ML 10,000bp")

print(xtable(mig750[,c(3:5,8:10)]), include.rownames=FALSE)



