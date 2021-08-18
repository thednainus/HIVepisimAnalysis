# function to sample from transmat
library(EpiModel)
library(lubridate)

sim <- readRDS("run1/results_sim.RDS")
sim_df <- as.data.frame(sim)

tm <- get_transmat(sim)

# total number of years simulated
years <- 40
#beginning of simulation time
init_sim_date <- ymd("1981-01-01")
# year of last sample dates in simulations
last_sample_date <- as.Date(x = years*365, origin = init_sim_date)



#read info from file and convert time from days to decimal years
art_init <- read.csv("run1/ART_init.csv")
art_init["time_decimal"] <- days2years(sampleTimes = art_init$time,
                                       init_date = init_sim_date)
dep <- read.csv("run1/departure_IDs.csv")
dep["time_decimal"] <- days2years(sampleTimes = dep$time,
                                  init_date = init_sim_date)

diag_info <- read.csv("run1/diag_time.csv")
diag_info["time_decimal"] <- days2years(sampleTimes = diag_info$time,
                                        init_date = init_sim_date)
stages <- read.csv("run1/stages.csv")
stages["time_decimal"] <- days2years(sampleTimes = stages$time,
                                     init_date = init_sim_date)



#sampled times
start_date <- ymd("1996-01-01")
start_date_dec <- decimal_date(start_date)
end_date <- ymd("2015-06-30")
end_date_dec <- decimal_date(end_date)

#total number of individuals in the transmission matrix
n_region <- get_n_by_origin(tm, location = "region")
n_global <- get_n_by_origin(tm, location = "global")

#get 5% of total number of individuals for each population
#nst = number of sample times to take for each population
nst_region <- round(0.05 * n_region,0)
nst_global <- round(0.05 * n_global,0)

#generate random date numbers from the start and end dates
st_ids_region <- sampleIDs(n_sample_times = nst_region, start_date = start_date_dec,
                           end_date = end_date_dec, art_init = art_init,
                           departure = dep, diag_info = diag_info,
                           tm = tm, location = "region")


st_ids_global <- sampleIDs(n_sample_times = nst_global,
                           start_date = decimal_date(init_sim_date),
                           end_date = decimal_date(last_sample_date),
                           art_init = art_init, departure = dep,
                           diag_info = diag_info,
                           tm = tm, location = "global")


rm(art_init, dep, diag_info, stages, start_date, start_date_dec, end_date, end_date_dec)
