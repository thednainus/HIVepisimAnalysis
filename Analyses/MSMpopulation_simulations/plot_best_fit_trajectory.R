#plot best fit trajectory
library(EpiModel)
library(HIVepisimAnalysis)
library(lubridate)
library(stringr)

dir_list <- dir("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/FitTrajectories/", full.names = T)

#beginning of simulation time
init_sim_date <- ymd("1980-01-01")

results <- data.frame()


for(i in 1:length(dir_list)){

  filename <- list.files(list.files(list.files
                                     (dir_list[i], full.names = T), full.names = T),
                          full.names = T, pattern = "results_sim.RDS")

  sim <- readRDS(filename)

  sim_param <- str_split(filename, "/")

  sim_df <- as.data.frame(sim)
  sim_df$param <- sim_param[[1]][11]
  sim_df$param <- unlist(lapply(sim_df$param, function(x) str_replace(string = x,
                                                           pattern = "small_msm_pop_sim",
                                                           replacement = "param")))
  sim_df$replicate <- sim_param[[1]][12]

  sim_df["rep_param"] <- paste(sim_df$replicate, sim_df$param, sep = "_")

  pop1_data_freq <- sim_df[c("rep_param", "time", "newDx_pop1")]
  pop1_data_freq["time_years"] <- days2years(pop1_data_freq$time, init_date = init_sim_date)

  pop1_data <- pop1_data_freq
  pop1_data["year"] <- unlist(lapply(pop1_data$time_years, function(x) strsplit(as.character(x), split = ".", fixed = TRUE)[[1]][1]))
  pop1_data$year <- as.factor(pop1_data$year)
  pop1_data$type <- as.factor(pop1_data$rep_param)

  pop1_data <- pop1_data[c(5, 1,3)]
  pop1_data_agg <- aggregate(newDx_pop1 ~ year + rep_param, data = pop1_data, FUN=sum)


  results <- rbind(results, pop1_data_agg)

}


#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))

#plot
results$year <- as.character(  results$year )
results$year <- as.numeric(  results$year )
sims = split( results , results$rep_param )
Y = do.call( cbind, lapply( sims , '[[', 'newDx_pop1' )  )
matplot( sims[[1]]$year, Y , type = 'l', col = 'grey', lwd=.5, ylim = c( 0, max(na.omit( c(results$newDx_pop1,  incidenceDiag$frequency)) ) ) )
#lines( sims[[ resample[1] ]]$year , sims[[ resample[1] ]]$newDx_pop1, col = 'red', lwd = 2 )
with( incidenceDiag, points( time, frequency, pty = 20, col = 'blue', cex = 1 ) )

