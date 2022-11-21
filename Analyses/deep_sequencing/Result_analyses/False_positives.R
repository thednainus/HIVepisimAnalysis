#Plot results for pairs identified as false positives with phyloscanner

library(stringr)
library(dplyr)

#read pairs in which filtering was W >= 80%

FP_W80 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_FP_results_W80_final.RDS")
FP_W80["param_mig"] <- paste(FP_W80$param, FP_W80$mig, sep = "_")
FP_W80["byGroup"] <- c(1:nrow(FP_W80))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)

total_chain <- FP_W80 %>% group_by(byGroup) %>%
  mutate(total_chain = length(str_split(chain, pattern = "-")[[1]]) - 2)

#get summary from total chain
summary_by_params <- total_chain %>%
  group_by(param_mig) %>%
  mutate(mean = mean(total_chain),
         median = median(total_chain),
         total_rows = n(),
         total_1_intermediate = sum(total_chain == 1)) %>%
  select(param_mig, mean, median, total_rows, total_1_intermediate) %>%
  distinct()


#read pairs in which filtering was W >= 80%

FP_W0.01 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_FP_results_W0.01_final.RDS")
FP_W0.01["param_mig"] <- paste(FP_W0.01$param, FP_W0.01$mig, sep = "_")
FP_W0.01["byGroup"] <- c(1:nrow(FP_W0.01))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)

total_chain_W0.01 <- FP_W0.01 %>% group_by(byGroup) %>%
  mutate(total_chain = length(str_split(chain, pattern = "-")[[1]]) - 2)

#get summary from total chain
summary_by_params_W0.01 <- total_chain_W0.01 %>%
  group_by(param_mig) %>%
  mutate(mean = mean(total_chain),
         median = median(total_chain),
         total_rows = n(),
         total_1_intermediate = sum(total_chain == 1)) %>%
  select(param_mig, mean, median, total_rows, total_1_intermediate) %>%
  distinct()
