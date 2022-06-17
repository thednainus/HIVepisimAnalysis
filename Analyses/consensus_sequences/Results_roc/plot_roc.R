#plot results from roc curves

#true trees
true_tree_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_treeTrees_s1.RDS")
true_tree_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_treeTrees_s2.RDS")

#ML tree: 1000bp
ml1000bp_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_ML1000bp_s1.RDS")
ml1000bp_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_ML1000bp_s2.RDS")

#ML tree: 10000bp
ml10000bp_s1 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_ML10000bp_s1.RDS")
ml10000bp_s2 <- readRDS("~/Box Sync/HIV_SanDiego/data_simulations/ROC_data_for_ploting/roc_ML10000bp_s2.RDS")

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

sampler1_1067_trueTrees <- subset(sampler1_1067, code == "TrueTrees")
sampler1_1067_1000bp <- subset(sampler1_1067, code == "1000bp")
sampler1_1067_10000bp <- subset(sampler1_1067, code == "10000bp")


#sampler 2
sampler2 <- rbind(true_tree_s2, ml1000bp_s2, ml10000bp_s2)
sampler2["param_mig_code"] <- paste(sampler2$param,
                                    sampler2$mig,
                                    sampler2$code,
                                    sep = "_")
sampler2_1067 <- subset(sampler2, param == "1067")
sampler2_1067 <- subset(sampler2, param == "2348")

#plot
quartz()
sampler1_1067 %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=param_mig_code)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))

quartz()
sampler1_1067_250mig %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=code)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_250migrants")

quartz()
sampler1_1067_500mig %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=code)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_500migrants")


quartz()
sampler1_1067_750mig %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=code)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_750migrants")

quartz()
sampler1_1067_trueTrees %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=mig)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_trueTrees")

quartz()
sampler1_1067_1000bp %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=mig)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_1000bp")


quartz()
sampler1_1067_10000bp %>%
  group_by(code, sampler, param, perc, mig, param_mig, specificity) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept= 0, linetype = 4) +
  geom_step(aes(color=mig)) +
  facet_wrap(~ perc, scales = "free", ncol = 4) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2)) +
  ggtitle("param1067_10000bp")


sampler1["param_mig_code_perc"] <- paste(sampler1$param,
                                         sampler1$mig,
                                         sampler1$code,
                                         sampler1$perc,
                                         sep = "_")

sampler1_AUC <- sampler1
sampler1_AUC["FPR"] <- 1 - sampler1_AUC$specificity
sampler1_AUC["TPR"] <- sampler1_AUC$sensitivity
auc <- sampler1_AUC %>%
  group_by(param, param_mig_code, mig, code, perc, param_mig_code_perc) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })

auc["mig_code"] <- paste(auc$mig,
                         auc$code,
                         sep = "_")

auc %>%
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
