#plot results from prc curves
library(DescTools)
library(caret)
library(dplyr)

prc_deepseq <- readRDS("prc_all_deepseq_w0.01.RDS")
prc_deepseq["param_mig"] <- paste(prc_deepseq$param, prc_deepseq$mig,
                                  sep = "_")

prc_deepseq_500 <- subset(prc_deepseq, mig == "500migrants")


prc_deepseq_250_750 <- subset(prc_deepseq, mig == "250migrants" |
                                mig == "750migrants")



#plot sampler 1, param 1067, 500mig ----
quartz()
prc_deepseq_500 %>%
  group_by(param, mig, param_mig) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = positives/(positives + negatives),
                 color = param), linetype = 4) +
  geom_line(aes(color = param), size = 0.8) +
  theme_bw() +
  ggtitle("PRC curves for results obtained with phyloscanner") +
  labs(y = 'Precision', x = 'Recall') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439'),
                      name = 'Parameters',
                      breaks = c('1067', '2348'),
                      labels = c('Combination 1', 'Combination 2'))




#Area under the curve ----


auc_W0.01 <- prc_deepseq_500 %>%
  group_by(param, mig, param_mig) %>%
  distinct() %>%
  group_modify(~ {
    AUC(.x$recall, .x$precision, na.rm = TRUE) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })



prc_deepseq_250_750$mig <- ifelse(prc_deepseq_250_750$mig == "250migrants",
                                  "mig = 1/0.68",
                                  "mig = 1/2.05")


quartz()
prc_deepseq_250_750 %>%
  group_by(param, mig, param_mig) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_hline(aes(yintercept = positives/(positives + negatives),
                 color = param), linetype = 4) +
  geom_line(aes(color = param), size = 0.8) +
  facet_wrap(~ mig) +
  theme_bw() +
  ggtitle("PRC curves for results obtained with phyloscanner") +
  labs(y = 'Precision', x = 'Recall') +
  theme(text = element_text(size = 20), legend.position = "bottom") +
  scale_colour_manual(values = c('#ad0075', '#9fc439'),
                      name = 'Parameters',
                      breaks = c('1067', '2348'),
                      labels = c('Combination 1', 'Combination 2'))




#Area under the curve ----


auc_W0.01_rest <- prc_deepseq_250_750 %>%
  group_by(param, mig, param_mig) %>%
  distinct() %>%
  group_modify(~ {
    AUC(.x$recall, .x$precision, na.rm = TRUE) %>%
      tibble::enframe(name = NULL, value = "AUC")
  })


