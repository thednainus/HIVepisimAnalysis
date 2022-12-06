# Summarize results for consensus sequences
# It will read the labels results (1 and 0 to determine if correctly detected
# a transmission pair or not)

library(tidyverse)
library(HIVepisimAnalysis)
library(caret)

get_prc_data <- function(df_data){

  thresholds <- seq(from = 0, to = 1, by = 0.01)
  rates <- lapply(thresholds, get_precision_recall, df_data)
  rates <- do.call(rbind, rates)

  return(rates)

}

# True Trees: sampler 1----
#for true trees
mig250_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_250mig.RDS")
mig250_true["mig"] <- "250"
mig500_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_500mig.RDS")
mig500_true["mig"] <- "500"
mig750_true <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_750mig.RDS")
mig750_true["mig"] <- "750"
true_trees <- rbind(mig250_true, mig500_true, mig750_true)

trueTrees_df <- data.frame(infectorProbability = true_trees$infectorProbability,
                           labels = true_trees$labels,
                           code = true_trees$Code,
                           sampler = true_trees$sampler,
                           param = true_trees$param,
                           perc = true_trees$perc,
                           mig = true_trees$mig,
                           rep = true_trees$rep,
                           param_mig = paste(true_trees$param,
                                             true_trees$mig,
                                             sep = "_"),
                           tag = paste(true_trees$Code,
                                       true_trees$sampler,
                                       true_trees$param,
                                       true_trees$perc,
                                       true_trees$mig,
                                       sep = "_"))

total_pairs_trueTrees <- trueTrees_df %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_mig, tag, total) %>%
  distinct()


prc_true_trees <- trueTrees_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })



saveRDS(prc_true_trees, "prc_trueTrees_s1.RDS")

rm(prc_true_trees)
# True Trees: sampler 2----
#for true trees
mig250_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_250mig.RDS")
mig250_true["mig"] <- "250"
mig500_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_500mig.RDS")
mig500_true["mig"] <- "500"
mig750_true <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_750mig.RDS")
mig750_true["mig"] <- "750"
true_trees <- rbind(mig250_true, mig500_true, mig750_true)



trueTrees_df <- data.frame(infectorProbability = true_trees$infectorProbability,
                           labels = true_trees$labels,
                           code = true_trees$Code,
                           sampler = true_trees$sampler,
                           param = true_trees$param,
                           perc = true_trees$perc,
                           mig = true_trees$mig,
                           rep = true_trees$rep,
                           param_mig = paste(true_trees$param,
                                             true_trees$mig,
                                             sep = "_"),
                           tag = paste(true_trees$Code,
                                       true_trees$sampler,
                                       true_trees$param,
                                       true_trees$perc,
                                       true_trees$mig,
                                       sep = "_"))

prc_true_trees <- trueTrees_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })

saveRDS(prc_true_trees, "prc_trueTrees_s2.RDS")



# ML trees 1000bp: sampler 1----
rm(prc_true_trees)

mig250 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_250mig.RDS")
mig250["mig"] <- "250"
mig500 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_500mig.RDS")
mig500["mig"] <- "500"
mig750 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_1000bp_750mig.RDS")
mig750["mig"] <- "750"
ml_1000bp <- rbind(mig250, mig500, mig750)




ml_1000bp_df <- data.frame(infectorProbability = ml_1000bp$infectorProbability,
                           labels = ml_1000bp$labels,
                           code = ml_1000bp$Code,
                           sampler = ml_1000bp$sampler,
                           param = ml_1000bp$param,
                           perc = ml_1000bp$perc,
                           mig = ml_1000bp$mig,
                           rep = ml_1000bp$rep,
                           param_mig = paste(ml_1000bp$param,
                                             ml_1000bp$mig,
                                             sep = "_"),
                           tag = paste(ml_1000bp$Code,
                                       ml_1000bp$sampler,
                                       ml_1000bp$param,
                                       ml_1000bp$perc,
                                       ml_1000bp$mig,
                                       sep = "_"))

total_pairs_1000bp <- ml_1000bp_df %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_mig, tag, total) %>%
  distinct()

prc_Ml1000bp <- ml_1000bp_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })



saveRDS(prc_Ml1000bp, "prc_ML1000bp_s1.RDS")


# ML trees 10,000bp: sampler 1----
rm(mig250, mig500, mig750, roc_Ml1000bp, ml_1000bp_df)
mig250 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_250mig.RDS")
mig250["mig"] <- "250"
mig500 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_500mig.RDS")
mig500["mig"] <- "500"
mig750 <- readRDS("Analyses/consensus_sequences/Results/sampler1/all_data_s1_10000bp_750mig.RDS")
mig750["mig"] <- "750"
ml_10000bp <- rbind(mig250, mig500, mig750)



ml_10000bp_df <- data.frame(infectorProbability = ml_10000bp$infectorProbability,
                            labels = ml_10000bp$labels,
                            code = ml_10000bp$Code,
                            sampler = ml_10000bp$sampler,
                            param = ml_10000bp$param,
                            perc = ml_10000bp$perc,
                            mig = ml_10000bp$mig,
                            rep = ml_10000bp$rep,
                            param_mig = paste(ml_10000bp$param,
                                              ml_10000bp$mig,
                                              sep = "_"),
                            tag = paste(ml_10000bp$Code,
                                        ml_10000bp$sampler,
                                        ml_10000bp$param,
                                        ml_10000bp$perc,
                                        ml_10000bp$mig,
                                        sep = "_"))

total_pairs_10000bp <- ml_10000bp_df %>% group_by(tag) %>%
  mutate(total = n()) %>%
  select(code, sampler, param, perc, mig, param_mig, tag, total) %>%
  distinct()

prc_Ml10000bp <- ml_10000bp_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })



saveRDS(prc_Ml10000bp, "prc_ML10000bp_s1.RDS")


# ML trees 1000bp: sampler 2----
rm(mig250, mig500, mig750, roc_Ml10000bp, ml_10000bp_df)
mig250 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_250mig.RDS")
mig250["mig"] <- "250"
mig500 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_500mig.RDS")
mig500["mig"] <- "500"
mig750 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_1000bp_750mig.RDS")
mig750["mig"] <- "750"
ml_1000bp <- rbind(mig250, mig500, mig750)




ml_1000bp_df <- data.frame(infectorProbability = ml_1000bp$infectorProbability,
                           labels = ml_1000bp$labels,
                           code = ml_1000bp$Code,
                           sampler = ml_1000bp$sampler,
                           param = ml_1000bp$param,
                           perc = ml_1000bp$perc,
                           mig = ml_1000bp$mig,
                           rep = ml_1000bp$rep,
                           param_mig = paste(ml_1000bp$param,
                                             ml_1000bp$mig,
                                             sep = "_"),
                           tag = paste(ml_1000bp$Code,
                                       ml_1000bp$sampler,
                                       ml_1000bp$param,
                                       ml_1000bp$perc,
                                       ml_1000bp$mig,
                                       sep = "_"))


prc_Ml1000bp <- ml_1000bp_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })



saveRDS(prc_Ml1000bp, "prc_ML1000bp_s2.RDS")


# ML trees 10000bp: sampler 2----
rm(mig250, mig500, mig750, roc_Ml1000bp, ml_1000bp_df)
mig250 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_250mig.RDS")
mig250["mig"] <- "250"
mig500 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_500mig.RDS")
mig500["mig"] <- "500"
mig750 <- readRDS("Analyses/consensus_sequences/Results/sampler2/all_data_s2_10000bp_750mig.RDS")
mig750["mig"] <- "750"
ml_10000bp <- rbind(mig250, mig500, mig750)



ml_10000bp_df <- data.frame(infectorProbability = ml_10000bp$infectorProbability,
                            labels = ml_10000bp$labels,
                            code = ml_10000bp$Code,
                            sampler = ml_10000bp$sampler,
                            param = ml_10000bp$param,
                            perc = ml_10000bp$perc,
                            mig = ml_10000bp$mig,
                            rep = ml_10000bp$rep,
                            param_mig = paste(ml_10000bp$param,
                                              ml_10000bp$mig,
                                              sep = "_"),
                            tag = paste(ml_10000bp$Code,
                                        ml_10000bp$sampler,
                                        ml_10000bp$param,
                                        ml_10000bp$perc,
                                        ml_10000bp$mig,
                                        sep = "_"))

prc_Ml10000bp <- ml_10000bp_df %>%
  group_by(tag) %>%
  group_map(~ {
    print(paste(.x$param[1], .x$perc[1], .x$rep[1], .x$mig[1], sep = "_"))
    #print(unique(.x$rep))
    get_prc_data(.x)

  })




saveRDS(prc_Ml10000bp, "prc_ML10000bp_s2.RDS")



