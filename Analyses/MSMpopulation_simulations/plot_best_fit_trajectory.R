#plot best fit trajectory
library(EpiModel)
library(HIVepisimAnalysis)
library(lubridate)
library(stringr)

#dir_list <- dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/FitTrajectories", full.names = T)
#dir_list <- dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/best_trajectories/MSMpop_25Oct2021", full.names = T)
#dir_list <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/results_v3", full.names = T)
#dir_list <- dir("/Users/user/Desktop/tmp2/osg_test_1000jobs/results_params_24Jan22", full.names = T)
dir_list <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results", full.names = T)
dir_list <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results_art_2004", full.names = T)
dir_list <- dir("/Users/user/Desktop/tmp2/osg_new_scrips/results_narrow_parameters", full.names = T)
#dir_list <- dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_50migrants", full.names = T)
#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

results_incid <- data.frame()
results_diag <- data.frame()

#params <- c("params_3158", "params_4299", "params_4013")
#params <- c("params_6582", "params_6876")
#params <- c("params_72", "params_9509")
#params <- c("params_5497", "params_2598")
#best params using the new scrips
#params <- c("params_113", "params_4352", "params_7007")
#params <- c("params_8598", "params_2575", "params_6284")
#art at 1995
params <- c("params_75", "params_308")
params <- c("params_2354", "params_887")
params <- c("params_3158", "params_887")
params <- c("params_5125", "params_887")

#art at 2005
params <- c("params_705", "params_802", "params_355")
params <- c("params_705", "params_802", "params_1352")
params <- c("params_1663", "params_1352")

#narrow parameter values (5000 jobs)
params <- c("params_64", "params_18")
params <- c("params_665", "params_118")
params <- c("params_665", "params_1067")

#"params_1067", "params_2348" was resampled for diagnosis
params <- c("params_665", "params_1067", "params_2348")

#params <- paste("/Users/user/Desktop/tmp2/osg_test_1000jobs/results_params_24Jan22", params, sep = "/")
params <- paste("/Users/user/Desktop/tmp2/osg_new_scrips/results", params, sep = "/")
params <- paste("/Users/user/Desktop/tmp2/osg_new_scrips/results_art_2004", params, sep = "/")
params <- paste("/Users/user/Desktop/tmp2/osg_new_scrips/results_narrow_parameters", params, sep = "/")

indexes <- unlist(lapply(params, function(x) which(dir_list == x)))
for(i in 1:length(indexes)){




  #filename <- list.files(list.files(list.files
  #                                   (dir_list[i], full.names = T), full.names = T),
  #                        full.names = T, pattern = "results_sim.RDS")
  index <- indexes[i]
  filename <- list.files(list.files(paste(dir_list[index], "rep_1", sep = "/"), full.names = T),
                          full.names = T, pattern = "results_sim.RDS")

  sim <- readRDS(filename)

  sim_param <- str_split(filename, "/")

  sim_df <- as.data.frame(sim)
  sim_df$param <- sim_param[[1]][8]
  #sim_df$param <- unlist(lapply(sim_df$param, function(x) str_replace(string = x,
  #                                                         replacement = "param")))
  #                                                         pattern = "small_msm_pop_sim",
  sim_df$replicate <- sim_param[[1]][9]

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


#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                 package = "HIVepisimAnalysis"))

#plot diagnosis
results_diag$year <- as.character(results_diag$year )
results_diag$year <- as.numeric(results_diag$year )
sims = split( results_diag , results_diag$rep_param )
Y = do.call( cbind, lapply( sims , '[[', 'newDx_pop1' )  )
matplot( sims[[1]]$year, Y[,1] , type = 'l', col = 'grey', lwd=.5,
         ylim = c( 0, max(na.omit( c(results_diag$newDx_pop1,  incidenceDiag$frequency)) ) ) )
#lines( sims[[ resample[1] ]]$year , sims[[ resample[1] ]]$newDx_pop1, col = 'red', lwd = 2 )
with( incidenceDiag, points( time, frequency, pty = 20, col = 'blue', cex = 1 ) )


#plot incidence
results_incid$year <- as.character(  results_incid$year )
results_incid$year <- as.numeric(  results_incid$year )
sims = split( results_incid , results_incid$rep_param )
Y = do.call( cbind, lapply( sims , '[[', 'incid.pop1' )  )
matplot( sims[[1]]$year, Y[,1] , type = 'l', col = 'grey', lwd=.5,
         ylim = c( 0, max(na.omit( c(results_incid$incid.pop1,  incidence$ECDC_incidence.N_Inf_M)) ) ) )
#lines( sims[[ resample[1] ]]$year , sims[[ resample[1] ]]$newDx_pop1, col = 'red', lwd = 2 )
with( incidence, points( ECDC_incidence.Year, ECDC_incidence.N_Inf_M, pty = 20, col = 'blue', cex = 1 ) )


library(ggplot2)
library(reshape2)

#merge data
#incidence
all_inc_data <- data.frame(year=incidence$ECDC_incidence.Year,
                           ECDC_incidence = incidence$ECDC_incidence.N_Inf_M,
                           best_fit_incidence1 = results_incid[1:41,3],
                           best_fit_incidence2 = results_incid[43:83, 3])

all_inc_data <- data.frame(year=incidence$ECDC_incidence.Year,
                           ECDC_incidence = incidence$ECDC_incidence.N_Inf_M,
                           best_fit_incidence1 = results_incid[1:41,3],
                           best_fit_incidence2 = results_incid[43:83, 3],
                           best_fit_incidence3 = results_incid[85:125, 3])

all_inc_data <- data.frame(year=incidence$ECDC_incidence.Year,
                           ECDC_incidence = incidence$ECDC_incidence.N_Inf_M,
                           best_fit_incidence1 = results_incid[1:41,3])


#melt data
incid <- melt(all_inc_data, id.vars = c("year"))

quartz()
ggplot(incid, aes(x = year, y = value)) +
  geom_line(aes(colour = variable), size = 1.5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  labs(y = "incidence")



#merge data
#incidence diagnosis
incidenceDiag$frequency
#all_dx_data <- data.frame(year = incidenceDiag$time,
#                          Diagnosis = incidenceDiag$frequency,
#                          best_fit_incid_diag = results_diag[43:84,3])

all_dx_data <- data.frame(year = incidenceDiag$time[6:41],
                          Diagnosis = incidenceDiag$frequency[6:41],
                          best_fit_incid_diag1 = results_diag[6:41,3],
                          best_fit_incid_diag2 = results_diag[48:83, 3])

all_dx_data <- data.frame(year = incidenceDiag$time[6:41],
                          Diagnosis = incidenceDiag$frequency[6:41],
                          best_fit_incid_diag1 = results_diag[6:41,3],
                          best_fit_incid_diag2 = results_diag[48:83, 3],
                          best_fit_incid_diag3 = results_diag[90:125, 3])

all_dx_data <- data.frame(year = incidenceDiag$time[6:41],
                          Diagnosis = incidenceDiag$frequency[6:41],
                          best_fit_incid_diag1 = results_diag[6:41,3],
                          best_fit_incid_diag2 = results_diag[48:83, 3],
                          best_fit_incid_diag1_2005 = results_2005[6:41,3],
                          best_fit_incid_diag2_2005 = results_2005[48:83, 3])


all_dx_data <- data.frame(year = incidenceDiag$time[6:41],
                          Diagnosis = incidenceDiag$frequency[6:41],
                          best_fit_incid_diag1 = results_diag[6:41,3])

best_fits <- data.frame(year = incidence$ECDC_incidence.Year,
                        best_fit_incid_diag1 = results_diag[1:41,3],
                        best_fit_incid_diag2 = results_diag[43:83,3],
                        best_fit_incidence1 = results_incid[1:41,3],
                        best_fit_incidence2 = results_incid[43:83,3]
                        )

best_fitsm <- melt(best_fits, id.vars = c("year"))

quartz()
ggplot(best_fitsm, aes(x = year, y = value)) +
  geom_line(aes(colour = variable), size = 1.5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20))


#melt data
diagn <- melt(all_dx_data, id.vars = c("year"))

quartz()
  ggplot(diagn, aes(x = year, y = value)) +
    geom_line(aes(colour = variable), size = 1.5) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 20)) +
    labs(y = "incid diagnosis")


params4dim <- readRDS(system.file("data/params4dim_14Jan2022.RDS",
                                  package = "HIVepisim"))[c(6582,6876),]

params4dim_all <- readRDS(system.file("data/params4dim_14Jan2022.RDS",
                                  package = "HIVepisim"))
#best incidence paramenter line number 6582
#best incidence diagnosis paramenter line number 6876


#plot different diagnosis
#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                 package = "HIVepisimAnalysis"))

diag_1995 <- readRDS("Results_for_diagnosis/diagnosis_1995_best887.RDS")
diag_1995_2 <- readRDS("Results_for_diagnosis/diagnosis_1995_best887diag_5125incid.RDS")
#diag_1995 <- diag_1995[48:83,]
diag_2004 <- readRDS("Results_for_diagnosis/diagnosis_2004.RDS")
diag_2004 <- diag_2004[48:83,]
diag_narrow <- readRDS("Results_for_diagnosis/diagnosis_narrow_parameters.RDS")
diag_narrow <- diag_narrow[48:83,]
diag_narrow2 <- readRDS("Results_for_diagnosis/diagnosis_narrow_parameters_best118.RDS")
diag_narrow2 <- diag_narrow2[48:83,]
diag_narrow3 <- readRDS("Results_for_diagnosis/diagnosis_narrow_params_best1067.RDS")
diag_narrow3 <- diag_narrow3[48:83,]

all_dx_data <- data.frame(year = incidenceDiag$time[6:41],
                          Diagnosis = incidenceDiag$frequency[6:41],
                          diag_1995_308 = diag_1995[c(6:41),3],
                          diag_1995_887 = diag_1995[c(48:83),3],
                          diag_1995_5125 = diag_1995_2[c(6:41), 3],
                          diag_narrow18 = diag_narrow[, 3],
                          diag_narrow118 = diag_narrow2[, 3],
                          diag_narrow1067 = diag_narrow3[, 3])


#melt data
diagn <- melt(all_dx_data, id.vars = c("year"))

quartz()
ggplot(diagn, aes(x = year, y = value)) +
  geom_line(aes(colour = variable), size = 1.5) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  labs(y = "incid diagnosis")




#nodes <- read.csv("total_nodes.csv")
#nodes <- rbind(c(1,30000,10000), nodes)
nodes1 <- rbind(c(1,60000,30000), nodes1)

params <- c("params_3025", "params_6552")
nodes1 <- read.csv("~/Desktop/tmp2/osg_new_scrips/results/params_3025/rep_1/results/total_nodes.csv")
nodes2 <- read.csv("~/Desktop/tmp2/osg_new_scrips/results/params_6552/rep_1/results/total_nodes.csv")

sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] , global_meandeg = (sim1$stats$nwstats$sim1[[1]][,2]*2)/nodes1$global)
sim1$stats$nwstats$sim1[[1]] <- cbind(sim1$stats$nwstats$sim1[[1]] ,
                                     region_meandeg = (sim$stats$nwstats$sim1[[1]][,4]*2)/nodes1$region)


plot(sim1, type = "formation", plots.joined = FALSE, stats = c("global_meandeg", "region_meandeg"))

plot(sim1, y=c("dall_pop1.flow", "dall_pop2.flow",
              "nArrivals_mig1", "nArrivals_mig2",
              "a1.flow", "a2.flow"), legend = TRUE)

plot(sim1, y=c("i.num.pop1", "s.num.pop1"), legend = TRUE)
plot(sim1, y=c("i.num.pop2", "s.num.pop2"), legend = TRUE)
plot(sim1, y=c("i.num.pop1", "s.num.pop1", "i.num.pop2", "s.num.pop2"), legend = TRUE)

plot(sim1, y=c("hstage0.pop1", "hstage1.pop1", "hstage2.pop1",
              "hstage3.pop1", "hstage.aids.pop1"), legend = TRUE)

plot(sim1, y=c("hstage0.pop2", "hstage1.pop2",
              "hstage2.pop2", "hstage3.pop2",
              "hstage.aids.pop2"), legend = TRUE)
