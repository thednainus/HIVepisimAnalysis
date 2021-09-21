#summarize epidata

library(DescTools)
library(lubridate)

#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

#epidata_dirs <- dir("simulation_results", full.names = T)
epidata_dirs <- dir("~/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/SanDiegoSimulationResults_smallpop", full.names = T)



epidata <- get_epi_data(epidata_dirs)

#get time, newDx and num.pop1 and i.num.pop1

pop1_data <- lapply(epidata, function(x) x[c("sim", "type", "time", "newDx_pop1", "num.pop1", "i.num.pop1")])
#pop1_data_perc <- lapply(pop1_data, function(x) data.frame(sim = x["sim"], type = x["type"],
#                                                           time = x["time"], newDx_perc = x["newDx_pop1"]/ x["num.pop1"]))

pop1_data_freq <- lapply(pop1_data, function(x) data.frame(sim = x["sim"], type = x["type"],
                                                           time = x["time"], newDx_pop1 = x["newDx_pop1"]))



pop1_data_freq <-do.call(rbind, pop1_data_freq)

#convert time from days to years
pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

pop1_data <- pop1_data_freq
pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
pop1_data$year <- as.factor(pop1_data$year)
pop1_data$type <- as.factor(pop1_data$type)
pop1_data$sim <- as.factor(pop1_data$sim)
pop1_data <- pop1_data[c(6, 1,2,4)]
pop1_data_agg <- aggregate(newDx_pop1 ~ year + sim + type, data = pop1_data, FUN=sum)
pop1_data_agg["param"] <- unlist(lapply(pop1_data_agg$type, function(x)
                                 str_replace(x, "small_msm_pop_sim", "param")))

#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
