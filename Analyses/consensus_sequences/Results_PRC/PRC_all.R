#plot results from PRC (precision recall curve)
library(DescTools)
library(caret)
library(DescTools)

#true trees
true_tree_s1 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_trueTrees_s1.RDS")
true_tree_s1 <- do.call(rbind, true_tree_s1)
true_tree_s2 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_trueTrees_s2.RDS")
true_tree_s2 <- do.call(rbind, true_tree_s2)

#ML tree: 1000bp
ml1000bp_s1 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_ML1000bp_s1.RDS")
ml1000bp_s1 <- do.call(rbind, ml1000bp_s1)

ml1000bp_s2 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_ML1000bp_s2.RDS")
ml1000bp_s2 <- do.call(rbind, ml1000bp_s2)

#ML tree: 10000bp
ml10000bp_s1 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_ML10000bp_s1.RDS")
ml10000bp_s1 <- do.call(rbind, ml10000bp_s1)

ml10000bp_s2 <- readRDS("Analyses/consensus_sequences/Results_PRC/results/prc_ML10000bp_s2.RDS")
ml10000bp_s2 <- do.call(rbind, ml10000bp_s2)


#sampler 1
sampler1 <- rbind(true_tree_s1, ml1000bp_s1, ml10000bp_s1)
sampler1["param_mig_code"] <- paste(sampler1$param,
                                    sampler1$mig,
                                    sampler1$code,
                                    sep = "_")
sampler1_1067 <- subset(sampler1, param == "1067")
sampler1_2348 <- subset(sampler1, param == "2348")



sampler1_1067_250mig <- subset(sampler1_1067, mig == "250")
sampler1_1067_500mig <- subset(sampler1_1067, mig == "500")
sampler1_1067_750mig <- subset(sampler1_1067, mig == "750")

sampler1_2348_250mig <- subset(sampler1_2348, mig == "250")
sampler1_2348_500mig <- subset(sampler1_2348, mig == "500")
sampler1_2348_750mig <- subset(sampler1_2348, mig == "750")




#sampler 2
sampler2 <- rbind(true_tree_s2, ml1000bp_s2, ml10000bp_s2)
sampler2["param_mig_code"] <- paste(sampler2$param,
                                    sampler2$mig,
                                    sampler2$code,
                                    sep = "_")
sampler2_1067 <- subset(sampler2, param == "1067")
sampler2_2348 <- subset(sampler2, param == "2348")

sampler2_1067_250mig <- subset(sampler2_1067, mig == "250")
sampler2_1067_500mig <- subset(sampler2_1067, mig == "500")
sampler2_1067_750mig <- subset(sampler2_1067, mig == "750")

sampler2_2348_250mig <- subset(sampler2_2348, mig == "250")
sampler2_2348_500mig <- subset(sampler2_2348, mig == "500")
sampler2_2348_750mig <- subset(sampler2_2348, mig == "750")




#Area under the curve ----

sampler1["param_mig_code_perc"] <- paste(sampler1$param,
                                         sampler1$mig,
                                         sampler1$code,
                                         sampler1$perc,
                                         sep = "_")


auc_sampler1 <- sampler1 %>%
  group_by(param, param_mig_code, mig, code, perc, param_mig_code_perc) %>%
  group_modify(~ {
    AUC(.x$recall, .x$precision, na.rm = TRUE) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })


sampler2["param_mig_code_perc"] <- paste(sampler2$param,
                                         sampler2$mig,
                                         sampler2$code,
                                         sampler2$perc,
                                         sep = "_")


auc_sampler2 <- sampler2 %>%
  group_by(param, param_mig_code, mig, code, perc, param_mig_code_perc) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })

auc_sampler2["mig_code"] <- paste(auc_sampler2$mig,
                                  auc_sampler2$code,
                                  sep = "_")




teste1 <- subset(auc_sampler1, param == "1067" & mig == "500")
teste1$sampler <- "1"
teste2 <- subset(auc_sampler2, param == "1067" & mig == "500")
teste2$sampler <- "2"

teste12 <- rbind(teste1, teste2)

teste12$mig_code <- as.factor(teste12$mig_code)
ggplot(data = teste12, aes(x = code, y = AUC)) +
  geom_bar(stat="identity", aes(col = code, linetype = sampler), fill = "white",
           size = 1.5,
           position=position_dodge()) +
  facet_wrap(~ perc, ncol = 5) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('TrueTrees', '1000bp',
                                 '10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp')) +
  coord_flip()


  group_by(mig_code, code) %>%


auc_sampler1 %>%
  group_by(mig_code, code) %>%
  ggplot(aes(x = mig_code, y = AUC)) +
  geom_bar(stat="count", aes(color = code, fill=code)) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  coord_flip()

auc_df <- as.data.frame(auc)
auc_df$AUC <- as.character(auc_df$AUC)
auc_df$AUC <- as.factor(auc_df$AUC)


ggplot(auc_df, aes(x = mig_code, y = AUC)) +
  geom_bar(stat="identity", aes(color = code, fill=code)) +
  facet_wrap(mig ~ perc, ncol = 5) +
  theme_bw()




teste <- subset(sampler1_AUC, param_mig_code_perc == "2348_250_10000bp_0.9")
