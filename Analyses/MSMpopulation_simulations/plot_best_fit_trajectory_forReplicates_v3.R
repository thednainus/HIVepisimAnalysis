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

#500 migrants
#dir_list <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/params72", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_500migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_500migrants/params_2348", full.names = T)


#50 migrants
#dir_list <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/params72", full.names = T)
dir_list <- dir("../Results_paper/older/best_trajectories_50migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/older/best_trajectories_50migrants/params_2348", full.names = T)

#75 migrants
dir_list <- dir("../Results_paper/best_trajectories_75migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_75migrants/params_2348", full.names = T)


#25 migrants
dir_list <- dir("../Results_paper/best_trajectories_25migrants/params_1067", full.names = T)
dir_list <- dir("../Results_paper/best_trajectories_250migrants/params_2348", full.names = T)

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

saveRDS(all_diag, "all_diag_m_and_q_2348_500migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_500migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_500migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_500migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_meandegree.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_75migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_75migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_1067_25migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_1067_25migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_meandegree.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_meandegree.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_75migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_75migrants.RDS")

saveRDS(all_diag, "all_diag_m_and_q_2348_25migrants.RDS")
saveRDS(all_incid, "all_incid_m_and_q_2348_25migrants.RDS")



#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                 package = "HIVepisimAnalysis"))

library(ggplot2)
library(reshape2)

#plot only the data in which migration rate was set to 50 migrants per year


param_1067 <- readRDS("Results_paper/all_diag_m_and_q_1067.RDS")
param_1067["param"] <- "1067"
param_1067["migrant"] <- "50"

param_2348 <- readRDS("Results_paper/all_diag_m_and_q_2348.RDS")
param_2348["param"] <- "2348"
param_2348["migrant"] <- "50"

param_1067_500 <- readRDS("Results_paper/all_diag_m_and_q_1067_500migrants.RDS")
param_1067_500["param"] <- "1067"
param_1067_500["migrant"] <- "500"

param_2348_500 <- readRDS("Results_paper/all_diag_m_and_q_2348_500migrants.RDS")
param_2348_500["param"] <- "2348"
param_2348_500["migrant"] <- "500"


all_diag <- rbind(param_1067[6:41,],
                  param_2348[6:41,],
                  param_1067_500[6:41,],
                  param_2348_500[6:41,])
all_diag["param_migrant"] <- paste(all_diag$param, all_diag$migrant, sep = "_")

all_diag$year <- as.character(all_diag$year)
all_diag$year <- as.numeric(all_diag$year)
all_diag$param_migrant <- as.factor(all_diag$param_migrant)

all_diag <- all_diag[,c(1:4,7)]

names(all_diag)[2:4] <- c("lower", "median", "upper")
all_diagm <- melt(all_diag, id.vars = c("year", "lower", "upper", "param_migrant"))


all_diag <- rbind(param_1067[6:41,],
                  param_2348[6:41,])
all_diag["param_migrant"] <- paste(all_diag$param, all_diag$migrant, sep = "_")

all_diag$year <- as.character(all_diag$year)
all_diag$year <- as.numeric(all_diag$year)
all_diag$param_migrant <- as.factor(all_diag$param_migrant)

all_diag <- all_diag[,c(1:4,7)]

names(all_diag)[2:4] <- c("lower", "median", "upper")
all_diagm <- melt(all_diag, id.vars = c("year", "lower", "upper", "param_migrant"))


quartz()
names(incidenceDiag)[1] <- "year"
incidenceDiag$year <- as.numeric(incidenceDiag$year)
ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  theme_bw() + ylab("Incidence of diagnosis") +
  scale_fill_manual(values=c("#c9222a", "#222ac9"), guide = "none") +
  scale_color_manual(values=c("#c9222a", "#222ac9", "black"),
                     breaks=c("1067_50", "2348_50", "San Diego data"),
                     labels=c("Parameters 1", "Parameters 2", "San Diego data")) +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(legend.title=element_blank())



ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  geom_line(data = incidenceDiag[6:41,], aes(y = frequency, colour = "San Diego data"),
            size = 0.8) +
  theme_bw() + ylab("Incidence of diagnosis") +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(legend.title=element_blank())


ggsave("incidence_diagnosis_50migrants.pdf", useDingbats=FALSE, width=12, height=9)



#Plot incidence of diagnoses for all the migration rates

param_1067_75 <- readRDS("Results_paper/all_diag_m_and_q_1067_75migrants.RDS")
param_1067_75["param"] <- "1067"
param_1067_75["migrant"] <- "75"

param_1067_25 <- readRDS("Results_paper/all_diag_m_and_q_1067_25migrants.RDS")
param_1067_25["param"] <- "1067"
param_1067_25["migrant"] <- "25"



param_2348_75 <- readRDS("Results_paper/all_diag_m_and_q_2348_75migrants.RDS")
param_2348_75["param"] <- "2348"
param_2348_75["migrant"] <- "75"


param_2348_25 <- readRDS("Results_paper/all_diag_m_and_q_2348_25migrants.RDS")
param_2348_25["param"] <- "2348"
param_2348_25["migrant"] <- "25"



all_diag <- rbind(param_1067[6:41,],
                  param_1067_25[6:41,],
                  param_1067_75[6:41,],
                  param_2348[6:41,],
                  param_2348_25[6:41,],
                  param_2348_75[6:41,])
all_diag["param_migrant"] <- paste(all_diag$param, all_diag$migrant, sep = "_")

all_diag$year <- as.character(all_diag$year)
all_diag$year <- as.numeric(all_diag$year)
all_diag$param_migrant <- as.factor(all_diag$param_migrant)

all_diag <- all_diag[,c(1:4,7)]


names(all_diag)[2:4] <- c("lower", "median", "upper")
all_diagm <- melt(all_diag, id.vars = c("year", "lower", "upper", "param_migrant"))


quartz()
names(incidenceDiag)[1] <- "year"
incidenceDiag$year <- as.numeric(incidenceDiag$year)
ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  geom_line(data = incidenceDiag[6:41,], aes(y = frequency, colour = "San Diego data"),
            size = 0.8) +
  theme_bw() + ylab("Incidence of diagnosis") +
  theme(text = element_text(size=20), legend.position = "bottom")





