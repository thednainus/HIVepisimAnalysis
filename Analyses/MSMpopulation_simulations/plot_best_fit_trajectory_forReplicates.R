#plot best fit trajectory
library(EpiModel)
library(HIVepisimAnalysis)
library(lubridate)
library(stringr)

m_and_qt <- function(dataframe){

  all_qt <- aggregate(dataframe[,3] ~ year, data = dataframe,
                      FUN=quantile, probs = c(0.025, 0.5, 0.975))

  all_qt <- data.frame(year = all_qt$year, all_qt$`dataframe[, 3]`)

  return(all_qt)


}


#dir_list <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/params72", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_50migrants/meandagree/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_50migrants/meandagree/params_2348", full.names = T)

#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

results_incid <- data.frame()
results_diag <- data.frame()

for(i in 1:length(dir_list)){



  #untar results files
  results_tar_files <- list.files(dir_list[i], full.names = T,
                                  pattern = "sim_results.tar.gz")
  untar(results_tar_files)

  filename <- "results/results_sim.RDS"

  sim <- readRDS(filename)

  sim_param <- str_split(results_tar_files, "/")

  sim_df <- as.data.frame(sim)
  sim_df$param <- str_split(sim_param[[1]][7], "/")[[1]][1]

  sim_df$replicate <- sim_param[[1]][8]

  sim_df["rep_param"] <- paste(sim_df$replicate, sim_df$param, sep = "_")

  #diagnosis
  pop1_data_freq <- sim_df[c("rep_param", "time", "newDx_pop1")]
  pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

  pop1_data <- pop1_data_freq
  pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
  pop1_data$year <- as.factor(pop1_data$year)
  pop1_data$type <- as.factor(pop1_data$rep_param)

  pop1_data <- pop1_data[c(5, 1,3)]
  pop1_data_agg <- aggregate(newDx_pop1 ~ year + rep_param, data = pop1_data, FUN=sum)

  #incidence
  pop1_data_incid <- sim_df[c("rep_param", "time", "incid.pop1")]
  pop1_data_incid["time_years"] <- days2years(pop1_data_incid$time, init_date = init_sim_date)

  pop1_data_incid2 <- pop1_data_incid
  pop1_data_incid2["year"] <- unlist(lapply(pop1_data_incid2$time_years,
                                            function(x) strsplit(as.character(x),
                                                                 split = ".", fixed = TRUE)[[1]][1]))
  pop1_data_incid2$year <- as.factor(pop1_data_incid2$year)
  pop1_data_incid2$type <- as.factor(pop1_data_incid2$rep_param)

  pop1_data_incid2 <- pop1_data_incid2[c(5, 1,3)]
  pop1_data_agg_incid <- aggregate(incid.pop1 ~ year + rep_param, data = pop1_data_incid2, FUN=sum)

  results_diag <- rbind(results_diag, pop1_data_agg)
  results_incid <- rbind(results_incid, pop1_data_agg_incid)

}

#get median and upper and lower quantiles
all_diag <- m_and_qt(results_diag)
all_incid <- m_and_qt(results_incid)

saveRDS(all_diag, "all_diag_m_and_q_1067_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_meandegree.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_meandegree.RDS")


#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                 package = "HIVepisimAnalysis"))

library(ggplot2)
library(reshape2)


param_1067 <- readRDS("all_diag_m_and_q_1067.RDS")
param_1067["param"] <- "1067"
param_1067_md <- readRDS("all_diag_m_and_q_1067_meandegree.RDS")
param_1067_md["param"] <- "1067_md"
param_2348 <- readRDS("all_diag_m_and_q_2348.RDS")
param_2348["param"] <- "2348"
param_2348_md <- readRDS("all_diag_m_and_q_2348_meandegree.RDS")
param_2348_md["param"] <- "2348_md"

#all_diag <- rbind(param_1067[6:41,], param_1067_md[6:41,],
#                  param_2348[6:41,], param_2348_md[6:41,])
all_diag <- rbind(param_1067[6:41,],
                  param_2348[6:41,])
#all_diag <- rbind(param_1067_md[6:41,],
#                  param_2348_md[6:41,])

params <- readRDS("Results_for_diagnosis/diagnosis_narrow_parameters_best1067_and2348.RDS")
params$year <- as.character(params$year)
params$year <- as.numeric(params$year)
params1067 <- params[48:83,]
params2348 <- params[90:125,]


names(all_diag)[2:4] <- c("lower", "median", "upper")
all_diagm <- melt(all_diag, id.vars = c("year", "lower", "upper", "param"))
all_diagm$year <- as.character(all_diagm$year)
all_diagm$year <- as.numeric(all_diagm$year)

quartz()
ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param), linetype = 2) +
  geom_line(data = incidenceDiag[6:41,], aes(y = frequency, colour = "San Diego data"),
            size = 0.8) +
  theme_bw() + ylab("Incidence of diagnosis") +
  scale_fill_manual(values=c("#c9222a", "#222ac9"), guide = "none") +
  scale_color_manual(values=c("#c9222a", "#222ac9", "black")) +
  theme(text = element_text(size=20), legend.position = "bottom")

names(incidenceDiag)[1] <- "year"



#incid_24Jan22 <- readRDS("params24Jan22_params82_incid.RDS")

#merge data
#incidence
#all_inc_data <- data.frame(year = all_incid$year[1:41],
#                           ECDC_incidence = incidence$ECDC_incidence.N_Inf_M,
#                           incid_24Jan22 = incid_24Jan22$incid.pop1[1:41],
#                           median = all_incid$X50.[1:41],
#                           upper = all_incid$X97.5.[1:41],
#                           lower = all_incid$X2.5.[1:41])

all_inc_data <- data.frame(year = all_incid$year[1:41],
                           ECDC_incidence = incidence$ECDC_incidence.N_Inf_M,
                           median = all_incid$X50.[1:41],
                           upper = all_incid$X97.5.[1:41],
                           lower = all_incid$X2.5.[1:41])


#melt data
incid <- melt(all_inc_data, id.vars = c("year", "lower", "upper"))
incid$year <- as.character(incid$year)
incid$year <- as.numeric(incid$year)


ggplot(incid, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.20) +
  geom_line(aes(y = value, linetype=variable)) +
  theme_bw()




#merge data
#incidence diagnosis
all_dx_data <- data.frame(year = param_1067$year[6:40],
                          Diagnosis = incidenceDiag$frequency[6:40],
                          median = all_diag$X50.[1:41],
                          upper = all_diag$X97.5.[1:41],
                          lower = all_diag$X2.5.[1:41])

#melt data
diagn <- melt(all_dx_data, id.vars = c("year", "upper", "lower"))
diagn$year <- as.character(diagn$year)
diagn$year <- as.numeric(diagn$year)

quartz()
ggplot(diagn, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.20) +
  geom_line(aes(y = value, linetype=variable)) +
  theme_bw()



params4dim <- readRDS(system.file("data/params4dim_14Jan2022.RDS",
                                  package = "HIVepisim"))[c(6582,6876),]

params4dim_all <- readRDS(system.file("data/params4dim_14Jan2022.RDS",
                                  package = "HIVepisim"))
#best incidence paramenter line number 6582
#best incidence diagnosis paramenter line number 6876