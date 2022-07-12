#plot results from roc curves
library(DescTools)
library(caret)


#true trees
true_tree_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_trueTrees_s1.RDS")
true_tree_s1 <- do.call(rbind, true_tree_s1)
true_tree_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_trueTrees_s2.RDS")
true_tree_s2 <- do.call(rbind, true_tree_s2)

#ML tree: 1000bp
ml1000bp_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_ML1000bp_s1.RDS")
ml1000bp_s1 <- do.call(rbind, ml1000bp_s1)

ml1000bp_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_ML1000bp_s2.RDS")
ml1000bp_s2 <- do.call(rbind, ml1000bp_s2)

#ML tree: 10000bp
ml10000bp_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_ML10000bp_s1.RDS")
ml10000bp_s1 <- do.call(rbind, ml10000bp_s1)

ml10000bp_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_myScript/merged_replicate_data/roc_ML10000bp_s2.RDS")
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


#plot sampler 1, param 1067, 500mig ----
quartz()
sampler1_1067_500mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 1 and 500 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_500_TrueTrees', '1067_500_1000bp',
                                 '1067_500_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 1, param 2348, 500mig ----
quartz()
sampler1_2348_500mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 1 and 500 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_500_TrueTrees', '2348_500_1000bp',
                                 '2348_500_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))



#plot sampler 1, param 1067, 250mig ----
quartz()
sampler1_1067_250mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 1 and 250 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_250_TrueTrees', '1067_250_1000bp',
                                 '1067_250_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 1, param 2348, 250mig ----
quartz()
sampler1_2348_250mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 1 and 250 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_250_TrueTrees', '2348_250_1000bp',
                                 '2348_250_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))


#plot sampler 1, param 1067, 750mig ----
quartz()
sampler1_1067_750mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 1 and 750 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_750_TrueTrees', '1067_750_1000bp',
                                 '1067_750_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 1, param 2348, 750mig ----
quartz()
sampler1_2348_750mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 1 and 750 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_750_TrueTrees', '2348_750_1000bp',
                                 '2348_750_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))




#plot sampler 2, param 1067, 500mig ----
quartz()
sampler2_1067_500mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 2 and 500 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_500_TrueTrees', '1067_500_1000bp',
                                 '1067_500_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 2, param 2348, 500mig ----
quartz()
sampler2_2348_500mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 2 and 500 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_500_TrueTrees', '2348_500_1000bp',
                                 '2348_500_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))



#plot sampler 2, param 1067, 250mig ----
quartz()
sampler2_1067_250mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 2 and 250 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_250_TrueTrees', '1067_250_1000bp',
                                 '1067_250_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 2, param 2348, 250mig ----
quartz()
sampler2_2348_250mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 2 and 250 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_250_TrueTrees', '2348_250_1000bp',
                                 '2348_250_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))


#plot sampler 2, param 1067, 750mig ----
quartz()
sampler2_1067_750mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 1 (sampler 2 and 750 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('1067_750_TrueTrees', '1067_750_1000bp',
                                 '1067_750_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))

#plot sampler 2, param 2348, 750mig ----
quartz()
sampler2_2348_750mig %>%
  group_by(code, sampler, param, perc, mig, param_mig_code) %>%
  ggplot(aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color = param_mig_code), size = 0.8) +
  facet_wrap(~ perc, scales = "free", ncol = 5) +
  theme_bw() +
  ggtitle("ROC curves for Parameters 2 (sampler 2 and 750 migrants)") +
  labs(y = 'True positive rate', x = 'False positive rate') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439', '#0187d3'),
                      name = 'Tree type',
                      breaks = c('2348_750_TrueTrees', '2348_750_1000bp',
                                 '2348_750_10000bp'),
                      labels = c('True trees', 'ML 1,000bp', 'ML 10,000bp'))


#Area under the curve ----

sampler1["param_mig_code_perc"] <- paste(sampler1$param,
                                         sampler1$mig,
                                         sampler1$code,
                                         sampler1$perc,
                                         sep = "_")

sampler1_AUC <- sampler1
sampler1_AUC["FPR"] <- 1 - sampler1_AUC$specificity
sampler1_AUC["TPR"] <- sampler1_AUC$sensitivity

auc_sampler1 <- sampler1_AUC %>%
  group_by(param, param_mig_code, mig, code, perc, param_mig_code_perc) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })

auc_sampler1["mig_code"] <- paste(auc_sampler1$mig,
                                  auc_sampler1$code,
                                  sep = "_")

sampler2["param_mig_code_perc"] <- paste(sampler2$param,
                                         sampler2$mig,
                                         sampler2$code,
                                         sampler2$perc,
                                         sep = "_")

sampler2_AUC <- sampler2
sampler2_AUC["FPR"] <- 1 - sampler2_AUC$specificity
sampler2_AUC["TPR"] <- sampler2_AUC$sensitivity

auc_sampler2 <- sampler2_AUC %>%
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
