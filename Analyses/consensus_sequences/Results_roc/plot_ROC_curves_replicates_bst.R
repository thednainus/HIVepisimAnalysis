#plot ROC curves


#sampler 1 ----

roc_true_trees_s1_50mig <- readRDS("all_true_roc_50migrants_sampler1_bst.RDS")
roc_true_trees_s1_50mig["param_perc"] <- paste(roc_true_trees_s1_50mig$param,
                                               roc_true_trees_s1_50mig$perc,
                                                sep = "_")

all_true_roc_s1_50mig_noNA <- na.omit(roc_true_trees_s1_50mig)

teste <- subset(all_true_roc_s1_50mig_noNA, rep == "1")


#roc_1000bp_s1 <- readRDS("all_1000bp_roc_sampler1_bst.RDS")
#roc_10000bp_s1 <- readRDS("all_10000bp_roc_sampler1_bst.RDS")

#all_data <- rbind(roc_true_trees_s1, roc_1000bp_s1, roc_10000bp_s1)
#all_data$code <- as.factor(all_data$code)
#all_data$code <- factor(all_data$code, levels = c("True trees", "1000", "10000"))

quartz()

teste_roc_0.05 <- subset(teste_roc, perc == "0.05")
teste_roc_0.9 <- subset(teste_roc, perc == "0.9")
quartz()
ggplot(data=teste_roc_0.9, aes(x=FPR, y=TPR, colour = param)) +
  geom_line() +
  geom_abline(slope= 1, intercept= 0, linetype = 4) + theme_bw(base_size = 20) +
  facet_wrap(~ perc, scales = "free", ncol = 4)+
  theme(legend.position="bottom") +
  ggtitle("Infector probability ROC curves: Diagnosed and not on ART") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate")



ggplot(data=roc_data, aes(x=FPR, y=TPR, colour = param_perc)) +
  geom_line(size = 1) +
  geom_abline(slope= 1, intercept= 0, linetype = 4) + theme_bw(base_size = 20) +
  facet_wrap(~ perc, scales = "free", ncol = 4)+
  theme(legend.position="bottom") +
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"),
                      name = "", labels = c("True trees",
                                            "Estimated trees: 1,000bp",
                                            "Estimated trees: 10,000bp")) +
  ggtitle("Infector probability ROC curves: Diagnosed and not on ART") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate")


#area under curve
auc <- all_data %>%
  group_by(perc, code) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })






#sampler 2 ----

roc_true_trees_s2 <- readRDS("all_roc_true_sampler2_bst.RDS")
roc_1000bp_s2 <- readRDS("all_1000bp_roc_s2_bst.RDS")
roc_10000bp_s2 <- readRDS("all_10000bp_roc_s2_bst.RDS")

all_data <- rbind(roc_true_trees_s2, roc_1000bp_s2, roc_10000bp_s2)
all_data$code <- as.factor(all_data$code)
all_data$code <- factor(all_data$code, levels = c("True trees", "1000", "10000"))


all_data$perc <- ifelse(all_data$perc == "5percSampler2", "5perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "10percSampler2", "10perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "20percSampler2", "20perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "30percSampler2", "30perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "40percSampler2", "40perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "50percSampler2", "50perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "60percSampler2", "60perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "70percSampler2", "70perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "80percSampler2", "80perc", all_data$perc)
all_data$perc <- ifelse(all_data$perc == "90percSampler2", "90perc", all_data$perc)

quartz()
ggplot(data=all_data, aes(x=FPR, y=TPR, colour = code)) +
  geom_line(size = 1) +
  geom_abline(slope= 1, intercept= 0, linetype = 4) + theme_bw(base_size = 20) +
  facet_wrap(~ perc, scales = "free", ncol = 4)+
  theme(legend.position="bottom") +
  scale_colour_manual(values=c("#0072B2", "#D55E00", "#009E73"),
                      name = "", labels = c("True trees",
                                            "Estimated trees: 1,000bp",
                                            "Estimated trees: 10,000bp")) +
  ggtitle("Infector probability ROC curves: Diagnosed and on ART or not on ART") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate")


#area under curve
auc <- all_data %>%
  group_by(perc, code) %>%
  group_modify(~ {
    AUC(.x$FPR, .x$TPR) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })
