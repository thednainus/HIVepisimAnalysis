#generate tables for paper

library(xtable)
source("Analyses/consensus_sequences/TP_TN_FP_FN_estimation.R")
source("Analyses/consensus_sequences/TP_TN_FP_FN_estimation_threshold90.R")

#Table for total pairs for consensus sequences

#for main text
#merge data for threshold of 80% and 90% for migration probability of 1/500
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_500mig_s1_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_500mig_s1_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s1_500mig_80_90 <- cbind(all_500mig_s1_80_data2[,c(2:4,8:13)],
                         all_500mig_s1_90_data2[,c(2:4,8:13)])

s1_500mig_80_90["TP"] <- paste(s1_500mig_80_90$total_TP80, s1_500mig_80_90$total_TP90, sep = "/")
s1_500mig_80_90["FP"] <- paste(s1_500mig_80_90$total_FP80, s1_500mig_80_90$total_FP90, sep = "/")
s1_500mig_80_90["TN"] <- paste(s1_500mig_80_90$total_TN80, s1_500mig_80_90$total_TN90, sep = "/")
s1_500mig_80_90["FN"] <- paste(s1_500mig_80_90$total_FN80, s1_500mig_80_90$total_FN90, sep = "/")
s1_500mig_80_90["Total"] <- paste(s1_500mig_80_90$total80, s1_500mig_80_90$total90, sep = "/")
s1_500mig_80_90["TP wiw"] <- paste(s1_500mig_80_90$total_wiw_within_TP80, s1_500mig_80_90$total_wiw_within_TP90, sep = "/")

final_table_500mig_s1 <- s1_500mig_80_90[c(1:6,37:42,55:60),c(1:3,19:22,24)]
final_table_500mig_s1 <- final_table_500mig_s1[order(final_table_500mig_s1$param,final_table_500mig_s1$perc),]
final_table_500mig_s1 <- final_table_500mig_s1[order(final_table_500mig_s1$code, decreasing = TRUE),]

final_table_500mig_s1["perc"] <- as.factor(as.numeric(final_table_500mig_s1$perc) * 100)

final_table_500mig_s1["param"] <- ifelse(final_table_500mig_s1$param == "1067",
                                         "Combination 1 ", "Combination 2 ")

final_table_500mig_s1["Parameter/Tree/Perc."] <- paste(final_table_500mig_s1$param,
                                                  final_table_500mig_s1$code,
                                                  final_table_500mig_s1$perc,
                                                  sep = "/ ")


print(xtable(final_table_500mig_s1[c(9,4:8)]), include.rownames=FALSE)



# sampler 1: for supplementary material ----
#merge data for threshold of 80% and 90% for migration probability of 1/250
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_250mig_s1_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_250mig_s1_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s1_250mig_80_90 <- cbind(all_250mig_s1_80_data2[,c(2:4,8:13)],
                         all_250mig_s1_90_data2[,c(2:4,8:13)])

s1_250mig_80_90["TP"] <- paste(s1_250mig_80_90$total_TP80, s1_250mig_80_90$total_TP90, sep = "/")
s1_250mig_80_90["FP"] <- paste(s1_250mig_80_90$total_FP80, s1_250mig_80_90$total_FP90, sep = "/")
s1_250mig_80_90["TN"] <- paste(s1_250mig_80_90$total_TN80, s1_250mig_80_90$total_TN90, sep = "/")
s1_250mig_80_90["FN"] <- paste(s1_250mig_80_90$total_FN80, s1_250mig_80_90$total_FN90, sep = "/")
s1_250mig_80_90["Total"] <- paste(s1_250mig_80_90$total80, s1_250mig_80_90$total90, sep = "/")
s1_250mig_80_90["TP wiw"] <- paste(s1_250mig_80_90$total_wiw_within_TP80, s1_250mig_80_90$total_wiw_within_TP90, sep = "/")

s1_250mig_80_90["perc"] <- as.factor(as.numeric(s1_250mig_80_90$perc) * 100)
s1_250mig_80_90["param"] <- ifelse(s1_250mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s1_250mig_80_90["Parameter/Tree/Perc."] <- paste(s1_250mig_80_90$param,
                                                 s1_250mig_80_90$code,
                                                 s1_250mig_80_90$perc,
                                                 sep = "/ ")


final_table_250mig_s1 <- s1_250mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_250mig_s1), include.rownames=FALSE)


#for supplementary material
#merge data for threshold of 80% and 90% for migration probability of 1/500
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_500mig_s1_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_500mig_s1_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s1_500mig_80_90 <- cbind(all_500mig_s1_80_data2[,c(2:4,8:13)],
                         all_500mig_s1_90_data2[,c(2:4,8:13)])

s1_500mig_80_90["TP"] <- paste(s1_500mig_80_90$total_TP80, s1_500mig_80_90$total_TP90, sep = "/")
s1_500mig_80_90["FP"] <- paste(s1_500mig_80_90$total_FP80, s1_500mig_80_90$total_FP90, sep = "/")
s1_500mig_80_90["TN"] <- paste(s1_500mig_80_90$total_TN80, s1_500mig_80_90$total_TN90, sep = "/")
s1_500mig_80_90["FN"] <- paste(s1_500mig_80_90$total_FN80, s1_500mig_80_90$total_FN90, sep = "/")
s1_500mig_80_90["Total"] <- paste(s1_500mig_80_90$total80, s1_500mig_80_90$total90, sep = "/")
s1_500mig_80_90["TP wiw"] <- paste(s1_500mig_80_90$total_wiw_within_TP80, s1_500mig_80_90$total_wiw_within_TP90, sep = "/")

s1_500mig_80_90["perc"] <- as.factor(as.numeric(s1_500mig_80_90$perc) * 100)
s1_500mig_80_90["param"] <- ifelse(s1_500mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s1_500mig_80_90["Parameter/Tree/Perc."] <- paste(s1_500mig_80_90$param,
                                                 s1_500mig_80_90$code,
                                                 s1_500mig_80_90$perc,
                                                 sep = "/ ")


final_table_500mig_s1 <- s1_500mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_500mig_s1), include.rownames=FALSE)


#for supplementary material
#merge data for threshold of 80% and 90% for migration probability of 1/750
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_750mig_s1_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_750mig_s1_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s1_750mig_80_90 <- cbind(all_750mig_s1_80_data2[,c(2:4,8:13)],
                         all_750mig_s1_90_data2[,c(2:4,8:13)])

s1_750mig_80_90["TP"] <- paste(s1_750mig_80_90$total_TP80, s1_750mig_80_90$total_TP90, sep = "/")
s1_750mig_80_90["FP"] <- paste(s1_750mig_80_90$total_FP80, s1_750mig_80_90$total_FP90, sep = "/")
s1_750mig_80_90["TN"] <- paste(s1_750mig_80_90$total_TN80, s1_750mig_80_90$total_TN90, sep = "/")
s1_750mig_80_90["FN"] <- paste(s1_750mig_80_90$total_FN80, s1_750mig_80_90$total_FN90, sep = "/")
s1_750mig_80_90["Total"] <- paste(s1_750mig_80_90$total80, s1_750mig_80_90$total90, sep = "/")
s1_750mig_80_90["TP wiw"] <- paste(s1_750mig_80_90$total_wiw_within_TP80, s1_750mig_80_90$total_wiw_within_TP90, sep = "/")

s1_750mig_80_90["perc"] <- as.factor(as.numeric(s1_750mig_80_90$perc) * 100)
s1_750mig_80_90["param"] <- ifelse(s1_750mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s1_750mig_80_90["Parameter/Tree/Perc."] <- paste(s1_750mig_80_90$param,
                                                 s1_750mig_80_90$code,
                                                 s1_750mig_80_90$perc,
                                                 sep = "/ ")


final_table_750mig_s1 <- s1_750mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_750mig_s1), include.rownames=FALSE)



# sampler 2: for supplementary material ----
#merge data for threshold of 80% and 90% for migration probability of 1/250
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_250mig_s2_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_250mig_s2_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s2_250mig_80_90 <- cbind(all_250mig_s2_80_data2[,c(2:4,8:13)],
                         all_250mig_s2_90_data2[,c(2:4,8:13)])

s2_250mig_80_90["TP"] <- paste(s2_250mig_80_90$total_TP80, s2_250mig_80_90$total_TP90, sep = "/")
s2_250mig_80_90["FP"] <- paste(s2_250mig_80_90$total_FP80, s2_250mig_80_90$total_FP90, sep = "/")
s2_250mig_80_90["TN"] <- paste(s2_250mig_80_90$total_TN80, s2_250mig_80_90$total_TN90, sep = "/")
s2_250mig_80_90["FN"] <- paste(s2_250mig_80_90$total_FN80, s2_250mig_80_90$total_FN90, sep = "/")
s2_250mig_80_90["Total"] <- paste(s2_250mig_80_90$total80, s2_250mig_80_90$total90, sep = "/")
s2_250mig_80_90["TP wiw"] <- paste(s2_250mig_80_90$total_wiw_within_TP80, s2_250mig_80_90$total_wiw_within_TP90, sep = "/")

s2_250mig_80_90["perc"] <- as.factor(as.numeric(s2_250mig_80_90$perc) * 100)
s2_250mig_80_90["param"] <- ifelse(s2_250mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s2_250mig_80_90["Parameter/Tree/Perc."] <- paste(s2_250mig_80_90$param,
                                                 s2_250mig_80_90$code,
                                                 s2_250mig_80_90$perc,
                                                 sep = "/ ")


final_table_250mig_s2 <- s2_250mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_250mig_s2), include.rownames=FALSE)


#for supplementary material
#merge data for threshold of 80% and 90% for migration probability of 1/500
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_500mig_s2_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_500mig_s2_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s2_500mig_80_90 <- cbind(all_500mig_s2_80_data2[,c(2:4,8:13)],
                         all_500mig_s2_90_data2[,c(2:4,8:13)])

s2_500mig_80_90["TP"] <- paste(s2_500mig_80_90$total_TP80, s2_500mig_80_90$total_TP90, sep = "/")
s2_500mig_80_90["FP"] <- paste(s2_500mig_80_90$total_FP80, s2_500mig_80_90$total_FP90, sep = "/")
s2_500mig_80_90["TN"] <- paste(s2_500mig_80_90$total_TN80, s2_500mig_80_90$total_TN90, sep = "/")
s2_500mig_80_90["FN"] <- paste(s2_500mig_80_90$total_FN80, s2_500mig_80_90$total_FN90, sep = "/")
s2_500mig_80_90["Total"] <- paste(s2_500mig_80_90$total80, s2_500mig_80_90$total90, sep = "/")
s2_500mig_80_90["TP wiw"] <- paste(s2_500mig_80_90$total_wiw_within_TP80, s2_500mig_80_90$total_wiw_within_TP90, sep = "/")

s2_500mig_80_90["perc"] <- as.factor(as.numeric(s2_500mig_80_90$perc) * 100)
s2_500mig_80_90["param"] <- ifelse(s2_500mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s2_500mig_80_90["Parameter/Tree/Perc."] <- paste(s2_500mig_80_90$param,
                                                 s2_500mig_80_90$code,
                                                 s2_500mig_80_90$perc,
                                                 sep = "/ ")


final_table_500mig_s2 <- s2_500mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_500mig_s2), include.rownames=FALSE)


#for supplementary material
#merge data for threshold of 80% and 90% for migration probability of 1/750
#individuals per year
#total pairs classified as TP, FP, TN, and FN
names(all_750mig_s2_80_data2)[c(5:13)] <- c("precision80",
                                            "sensitivity80",
                                            "specificity80",
                                            "total_TP80",
                                            "total_FP80",
                                            "total_TN80",
                                            "total_FN80",
                                            "total80",
                                            "total_wiw_within_TP80")

names(all_750mig_s2_90_data2)[c(5:13)] <- c("precision90",
                                            "sensitivity90",
                                            "specificity90",
                                            "total_TP90",
                                            "total_FP90",
                                            "total_TN90",
                                            "total_FN90",
                                            "total90",
                                            "total_wiw_within_TP90")

s2_750mig_80_90 <- cbind(all_750mig_s2_80_data2[,c(2:4,8:13)],
                         all_750mig_s2_90_data2[,c(2:4,8:13)])

s2_750mig_80_90["TP"] <- paste(s2_750mig_80_90$total_TP80, s2_750mig_80_90$total_TP90, sep = "/")
s2_750mig_80_90["FP"] <- paste(s2_750mig_80_90$total_FP80, s2_750mig_80_90$total_FP90, sep = "/")
s2_750mig_80_90["TN"] <- paste(s2_750mig_80_90$total_TN80, s2_750mig_80_90$total_TN90, sep = "/")
s2_750mig_80_90["FN"] <- paste(s2_750mig_80_90$total_FN80, s2_750mig_80_90$total_FN90, sep = "/")
s2_750mig_80_90["Total"] <- paste(s2_750mig_80_90$total80, s2_750mig_80_90$total90, sep = "/")
s2_750mig_80_90["TP wiw"] <- paste(s2_750mig_80_90$total_wiw_within_TP80, s2_750mig_80_90$total_wiw_within_TP90, sep = "/")

s2_750mig_80_90["perc"] <- as.factor(as.numeric(s2_750mig_80_90$perc) * 100)
s2_750mig_80_90["param"] <- ifelse(s2_750mig_80_90$param == "1067",
                                   "Combination 1 ", "Combination 2 ")

s2_750mig_80_90["Parameter/Tree/Perc."] <- paste(s2_750mig_80_90$param,
                                                 s2_750mig_80_90$code,
                                                 s2_750mig_80_90$perc,
                                                 sep = "/ ")


final_table_750mig_s2 <- s2_750mig_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                           45:43,51:49,57:55,6:4, 12:10,18:16,
                                           24:22,30:28,36:34,42:40,48:46,54:52,
                                           60:58),c(25,19:22,24)]



print(xtable(final_table_750mig_s2), include.rownames=FALSE)




#for sampler 1: supplementary material: sensitivity, specifictiy and precision


s1_250mig_pss_80_90 <- cbind(all_250mig_s1_80_data2[,c(2:4,5:7)],
                         all_250mig_s1_90_data2[,c(2:4,5:7)])

s1_250mig_pss_80_90["Sensitivity"] <- paste(round(s1_250mig_pss_80_90$sensitivity80, 3),
                                            round(s1_250mig_pss_80_90$sensitivity90, 3), sep = "/")
s1_250mig_pss_80_90["Specificity"] <- paste(round(s1_250mig_pss_80_90$specificity80, 3),
                                            round(s1_250mig_pss_80_90$specificity90, 3), sep = "/")
s1_250mig_pss_80_90["Precision"] <- paste(round(s1_250mig_pss_80_90$precision80, 3),
                                          round(s1_250mig_pss_80_90$precision90, 3), sep = "/")

final_table_250mig_pss_s1 <- s1_250mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_250mig_pss_s1["perc"] <- as.factor(as.numeric(final_table_250mig_pss_s1$perc) * 100)

final_table_250mig_pss_s1["param"] <- ifelse(final_table_250mig_pss_s1$param == "1067",
                                         "Combination 1 ", "Combination 2 ")

final_table_250mig_pss_s1["Parameter/Tree/Perc."] <- paste(final_table_250mig_pss_s1$param,
                                                           final_table_250mig_pss_s1$code,
                                                           final_table_250mig_pss_s1$perc,
                                                  sep = "/ ")


print(xtable(final_table_250mig_pss_s1[c(7,4:6)]), include.rownames=FALSE)


#500mig
s1_500mig_pss_80_90 <- cbind(all_500mig_s1_80_data2[,c(2:4,5:7)],
                             all_500mig_s1_90_data2[,c(2:4,5:7)])

s1_500mig_pss_80_90["Sensitivity"] <- paste(round(s1_500mig_pss_80_90$sensitivity80, 3),
                                            round(s1_500mig_pss_80_90$sensitivity90, 3), sep = "/")
s1_500mig_pss_80_90["Specificity"] <- paste(round(s1_500mig_pss_80_90$specificity80, 3),
                                            round(s1_500mig_pss_80_90$specificity90, 3), sep = "/")
s1_500mig_pss_80_90["Precision"] <- paste(round(s1_500mig_pss_80_90$precision80, 3),
                                          round(s1_500mig_pss_80_90$precision90, 3), sep = "/")

final_table_500mig_pss_s1 <- s1_500mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_500mig_pss_s1["perc"] <- as.factor(as.numeric(final_table_500mig_pss_s1$perc) * 100)

final_table_500mig_pss_s1["param"] <- ifelse(final_table_500mig_pss_s1$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_500mig_pss_s1["Parameter/Tree/Perc."] <- paste(final_table_500mig_pss_s1$param,
                                                           final_table_500mig_pss_s1$code,
                                                           final_table_500mig_pss_s1$perc,
                                                           sep = "/ ")


print(xtable(final_table_500mig_pss_s1[c(7,4:6)]), include.rownames=FALSE)


#750mig
s1_750mig_pss_80_90 <- cbind(all_750mig_s1_80_data2[,c(2:4,5:7)],
                             all_750mig_s1_90_data2[,c(2:4,5:7)])

s1_750mig_pss_80_90["Sensitivity"] <- paste(round(s1_750mig_pss_80_90$sensitivity80, 3),
                                            round(s1_750mig_pss_80_90$sensitivity90, 3), sep = "/")
s1_750mig_pss_80_90["Specificity"] <- paste(round(s1_750mig_pss_80_90$specificity80, 3),
                                            round(s1_750mig_pss_80_90$specificity90, 3), sep = "/")
s1_750mig_pss_80_90["Precision"] <- paste(round(s1_750mig_pss_80_90$precision80, 3),
                                          round(s1_750mig_pss_80_90$precision90, 3), sep = "/")

final_table_750mig_pss_s1 <- s1_750mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_750mig_pss_s1["perc"] <- as.factor(as.numeric(final_table_750mig_pss_s1$perc) * 100)

final_table_750mig_pss_s1["param"] <- ifelse(final_table_750mig_pss_s1$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_750mig_pss_s1["Parameter/Tree/Perc."] <- paste(final_table_750mig_pss_s1$param,
                                                           final_table_750mig_pss_s1$code,
                                                           final_table_750mig_pss_s1$perc,
                                                           sep = "/ ")


print(xtable(final_table_750mig_pss_s1[c(7,4:6)]), include.rownames=FALSE)


#for sampler 2: supplementary material: sensitivity, specifictiy and precision ----


s2_250mig_pss_80_90 <- cbind(all_250mig_s2_80_data2[,c(2:4,5:7)],
                             all_250mig_s2_90_data2[,c(2:4,5:7)])

s2_250mig_pss_80_90["Sensitivity"] <- paste(round(s2_250mig_pss_80_90$sensitivity80, 3),
                                            round(s2_250mig_pss_80_90$sensitivity90, 3), sep = "/")
s2_250mig_pss_80_90["Specificity"] <- paste(round(s2_250mig_pss_80_90$specificity80, 3),
                                            round(s2_250mig_pss_80_90$specificity90, 3), sep = "/")
s2_250mig_pss_80_90["Precision"] <- paste(round(s2_250mig_pss_80_90$precision80, 3),
                                          round(s2_250mig_pss_80_90$precision90, 3), sep = "/")

final_table_250mig_pss_s2 <- s2_250mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_250mig_pss_s2["perc"] <- as.factor(as.numeric(final_table_250mig_pss_s2$perc) * 100)

final_table_250mig_pss_s2["param"] <- ifelse(final_table_250mig_pss_s2$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_250mig_pss_s2["Parameter/Tree/Perc."] <- paste(final_table_250mig_pss_s2$param,
                                                           final_table_250mig_pss_s2$code,
                                                           final_table_250mig_pss_s2$perc,
                                                           sep = "/ ")


print(xtable(final_table_250mig_pss_s2[c(7,4:6)]), include.rownames=FALSE)


#500mig
s2_500mig_pss_80_90 <- cbind(all_500mig_s2_80_data2[,c(2:4,5:7)],
                             all_500mig_s2_90_data2[,c(2:4,5:7)])

s2_500mig_pss_80_90["Sensitivity"] <- paste(round(s2_500mig_pss_80_90$sensitivity80, 3),
                                            round(s2_500mig_pss_80_90$sensitivity90, 3), sep = "/")
s2_500mig_pss_80_90["Specificity"] <- paste(round(s2_500mig_pss_80_90$specificity80, 3),
                                            round(s2_500mig_pss_80_90$specificity90, 3), sep = "/")
s2_500mig_pss_80_90["Precision"] <- paste(round(s2_500mig_pss_80_90$precision80, 3),
                                          round(s2_500mig_pss_80_90$precision90, 3), sep = "/")

final_table_500mig_pss_s2 <- s2_500mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_500mig_pss_s2["perc"] <- as.factor(as.numeric(final_table_500mig_pss_s2$perc) * 100)

final_table_500mig_pss_s2["param"] <- ifelse(final_table_500mig_pss_s2$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_500mig_pss_s2["Parameter/Tree/Perc."] <- paste(final_table_500mig_pss_s2$param,
                                                           final_table_500mig_pss_s2$code,
                                                           final_table_500mig_pss_s2$perc,
                                                           sep = "/ ")


print(xtable(final_table_500mig_pss_s2[c(7,4:6)]), include.rownames=FALSE)


#750mig
s2_750mig_pss_80_90 <- cbind(all_750mig_s2_80_data2[,c(2:4,5:7)],
                             all_750mig_s2_90_data2[,c(2:4,5:7)])

s2_750mig_pss_80_90["Sensitivity"] <- paste(round(s2_750mig_pss_80_90$sensitivity80, 3),
                                            round(s2_750mig_pss_80_90$sensitivity90, 3), sep = "/")
s2_750mig_pss_80_90["Specificity"] <- paste(round(s2_750mig_pss_80_90$specificity80, 3),
                                            round(s2_750mig_pss_80_90$specificity90, 3), sep = "/")
s2_750mig_pss_80_90["Precision"] <- paste(round(s2_750mig_pss_80_90$precision80, 3),
                                          round(s2_750mig_pss_80_90$precision90, 3), sep = "/")

final_table_750mig_pss_s2 <- s2_750mig_pss_80_90[c(3:1,9:7,15:13,21:19,27:25,33:31,39:37,
                                                   45:43,51:49,57:55,6:4, 12:10,18:16,
                                                   24:22,30:28,36:34,42:40,48:46,54:52,
                                                   60:58),c(1:3,13:15)]

final_table_750mig_pss_s2["perc"] <- as.factor(as.numeric(final_table_750mig_pss_s2$perc) * 100)

final_table_750mig_pss_s2["param"] <- ifelse(final_table_750mig_pss_s2$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_750mig_pss_s2["Parameter/Tree/Perc."] <- paste(final_table_750mig_pss_s2$param,
                                                           final_table_750mig_pss_s2$code,
                                                           final_table_750mig_pss_s2$perc,
                                                           sep = "/ ")


print(xtable(final_table_750mig_pss_s2[c(7,4:6)]), include.rownames=FALSE)





#for main text
s1_250mig_pss_80_90 <- cbind(all_250mig_s1_80_data2[,c(2:4,5:7)],
                             all_250mig_s1_90_data2[,c(2:4,5:7)])

s1_250mig_pss_80_90["Sensitivity"] <- paste(round(s1_250mig_pss_80_90$sensitivity80, 3),
                                            round(s1_250mig_pss_80_90$sensitivity90, 3), sep = "/")
s1_250mig_pss_80_90["Specificity"] <- paste(round(s1_250mig_pss_80_90$specificity80, 3),
                                            round(s1_250mig_pss_80_90$specificity90, 3), sep = "/")
s1_250mig_pss_80_90["Precision"] <- paste(round(s1_250mig_pss_80_90$precision80, 3),
                                          round(s1_250mig_pss_80_90$precision90, 3), sep = "/")

final_table_250mig_pss_s1 <- s1_250mig_pss_80_90[c(3:1,39:37,57:55,6:4,42:40,60:58),c(1:3,13:15)]

final_table_250mig_pss_s1["perc"] <- as.factor(as.numeric(final_table_250mig_pss_s1$perc) * 100)

final_table_250mig_pss_s1["param"] <- ifelse(final_table_250mig_pss_s1$param == "1067",
                                             "Combination 1 ", "Combination 2 ")

final_table_250mig_pss_s1["Parameter/Tree/Perc."] <- paste(final_table_250mig_pss_s1$param,
                                                           final_table_250mig_pss_s1$code,
                                                           final_table_250mig_pss_s1$perc,
                                                           sep = "/ ")


print(xtable(final_table_250mig_pss_s1[c(7,4:6)]), include.rownames=FALSE)



#500mig




