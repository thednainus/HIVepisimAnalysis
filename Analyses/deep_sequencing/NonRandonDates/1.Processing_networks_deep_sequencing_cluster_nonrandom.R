# Get transmission matrix ----
# Script has been checked on 11 July 2022
# In this script I get 1 pair that represent a transmission pair in which
# sample date of susceptible individual is time of infection + 1 week
# and run VirusTreeSimulator

# This code has been tested on a Mac OS and might not work on windows or linux

library(EpiModel)
library(HIVepisim)
library(HIVepisimAnalysis)
library(DescTools)
library(stringr)
library(ape)
library(phydynR)
library(castor)
library(dplyr)
library(lubridate)

param_list <- commandArgs(trailingOnly = TRUE)

# percentage of population to sampled IDs
line_number <- as.numeric(param_list[1])

seed_value <- as.numeric(param_list[2])
message(seed_value)
#seed_value <- 589634

set.seed(seed_value)

params_all <- readRDS(system.file("data/simulations_imperial_cluster.RDS",
                                  package = "HIVepisimAnalysis"))

#params_1067
params_all <- params_all[c(901:1000),]

#params_2348
#params_all <- params_all[c(1901:1930),]

params <- params_all[line_number,]

perc_pop_region <- params$perc


#untar tar file
#simulation_data <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants/params_1067/rep_12/sim_results.tar.gz"
simulation_data <- paste(getwd(), "sim_results.tar.gz", sep = "/")
untar(simulation_data)


# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output_deepseq" if directory does not exist
# and will save results to this directory "output_deepseq"

# Location for VirusTreeSimulator. It should be changed to the correct location on your computer.
#Software <- "java -jar VirusTreeSimulator/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
#parameters <- "-demoModel Constant -N0 1"
# parameters following Ratman et al. 2016 (Mol Biol Evol 34: 185-203)
# parameters in Ratman et al. 2016 is per year, and below I converted it to
# days as my simulations are in units of days
# effective population per year growth rate (r) =  2.851904
#  -t50 The time point, relative to the time of infection in backwards time, at
# which the population is equal to half its final asymptotic value,
# t50 in Ratman et al paper = -2 years
parameters <- "-demoModel Logistic -N0 1 -growthRate 2.851904 -t50 -2 -forceCoalescence"

#total number of years simulated
years <-  41


# year of last sample dates in simulations
#beginning of simulation time
#first case of HIV in San Diego was reported in 1980
init_sim_date <- ymd("1980-01-01")
# end of simulation date
end_sim_date <- ymd("2021-01-01")
# year of last sample dates in simulations
last_sample_date <- end_sim_date


#Times for sampling IDs and sampling times
#sample individuals in the past 5 years
start_date <- ymd("2015-01-01")
start_date_dec <- decimal_date(start_date)
end_date_dec <- decimal_date(last_sample_date)


#Create directory named output_deepseq if it does not exist
if (!dir.exists("output_deepseq")) {
  dir.create("output_deepseq")
}

# read departure ID files
#read info from file and convert time from days to decimal years
art_init <- read.csv("results/ART_init.csv")
art_init["time_decimal"] <- days2years(sampleTimes = art_init$time,
                                       init_date = init_sim_date)
dep <- read.csv("results/departure_IDs.csv")
dep["time_decimal"] <- days2years(sampleTimes = dep$time,
                                  init_date = init_sim_date)
dep["IDs"] <- dep$infID

diag_info <- read.csv("results/diag_time.csv")
diag_info["time_decimal"] <- days2years(sampleTimes = diag_info$time,
                                        init_date = init_sim_date)

stages <- read.csv("results/stages.csv")
stages["time_decimal"] <- days2years(sampleTimes = stages$time,
                                     init_date = init_sim_date)


origin <- read.csv("results/origin.csv")
origin["time_decimal"] <- days2years(sampleTimes = origin$time,
                                     init_date = init_sim_date)

#read simulation results
sim <- readRDS("results/results_sim.RDS")
sim_df <- as.data.frame(sim)

tm <- get_transmat(sim)
tm["year"] <-days2years(tm$at, init_date = init_sim_date)


# create file to be used with VirusTreeSimulator
output <- paste("output_deepseq", "output", sep ="/")
#output <- paste(output, PLWHIV, newinf, sep = "_")

#inf file name for VirusTreeSimulator
inf_file <- paste(output, "_inf.csv", sep = "")
#sample file name for VirusTreeSimulator
sample_file <- paste(output, "_sample.csv", sep = "")

seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

create_inf_csv(tm, time_tr = rep(1980, length(seed_names)), prefix = output)

# sample IDs and time of sampling
st_ids_region_all <- sampleIDs(perc = 0.1,
                               start_date = start_date_dec,
                               end_date = end_date_dec,
                               art_init = art_init,
                               departure = dep,
                               diag_info = diag_info,
                               origin = origin,
                               tm = tm,
                               location = "region")

#write proportion of sampled ids that were not in region anymore
prop_not_region <-
  nrow(st_ids_region_all[(st_ids_region_all$migrant == 2 |
                            st_ids_region_all$migrant == 12), ]) / nrow(st_ids_region_all)

#remove rows in which individuals are not in region anymore
#before or at time of sampling
st_ids_region <-  st_ids_region_all[(st_ids_region_all$migrant == 1 |
                       st_ids_region_all$migrant == 21), ]


list_samplesToAnalyse <- data.frame()

#get 10 pairs the sample time of sus is the time of infection
#then get just the first pair
for (i in 1:15) {
  print(i)
  tm_sampled <- tm[tm$inf == st_ids_region$sampled_ID[i], ]
  if (nrow(tm_sampled) != 0) {
    pair_sus <- sample_n(tm_sampled, 1)
    migrant_sus <- get_origin_at_samplingTime(tm, origin,
                                              pair_sus$sus,
                                              (pair_sus$year + 0.01))

    if((migrant_sus == 1 | migrant_sus == 21)){
      new_data <-
        data.frame(
          sampled_ID = c(st_ids_region$sampled_ID[i], pair_sus$sus),
          sampled_time = c(st_ids_region$sampled_time[i], (pair_sus$year + 0.01)),
          migrant = c(st_ids_region$migrant[i], migrant_sus)
        )

      list_samplesToAnalyse <- rbind(list_samplesToAnalyse, new_data)
    }

  }
}

# sampled IDs are not on ART
# Sample 10 viruses per ID
# note that this is just to start the simulations
# but I should try to find real values for those

#to select just one pair to analyse with phyloscanner
st_ids_region <- list_samplesToAnalyse[1:2,]

#st_ids_region["date"] <-
#  as.Date(date_decimal(st_ids_region$sampled_time))
#st_ids_region["time_days"] <-
#  as.numeric(st_ids_region$date - init_sim_date)
#st_ids_region["time_years"] <-
#  st_ids_region$time_days * (1 / 365)
st_ids_region["tip_name"] <- paste(st_ids_region$sampled_ID,
                                   st_ids_region$migrant,
                                    sep = "_")

create_sample_csv2(
  ids = st_ids_region$sampled_ID,
  time_seq = st_ids_region$sampled_time,
  seq_count = 11,
  prefix = output
)

#Create directory named VTS (for VirusTreeSimulator) if it does not exist
if (!dir.exists("output_deepseq/vts/")) {
  dir.create("output_deepseq/vts/")
}

#prefix with location of output directory to save results of VirusTreeSimulator
prefix_vts <-  "output_deepseq/vts/results_vts"

# Run VirusTreeSimulator
cmd <- paste(Software,
             parameters,
             inf_file,
             sample_file,
             prefix_vts, sep = " ")
system(cmd)

#Create directory to save phylogenetic trees
if (!dir.exists("output_deepseq/vts/merged_trees/")) {
  dir.create("output_deepseq/vts/merged_trees")
}

#prefix with location of output directory to save merged phylogenetic trees
prefix_trees <-  "output_deepseq/vts/merged_trees/results"

#read all trees from vts
list_trees <- dir("output_deepseq/vts", pattern = "*_detailed.nex", full.names = TRUE)
trees <- lapply(list_trees, read.nexus)
#add root.edge
if(length(trees) > 1){
  trees_rootedge <- lapply(trees, add_root_edge2,
                           root.edge_value = 0)
  vts_tree <- merge_trees(trees_rootedge)
} else{
  trees_rootedge <- add_root_edge2(tree = trees[[1]], root.edge_value = 0)
  vts_tree <- trees[[1]]
  }



#now remove only sampled_11. I simulated 11 sequences per ID
#because the detailed version of the virustreesimulator
#output file when read into R generated a weird number of internal numbers
#but removing some sequences, it fixes the number of internal nodes in the tree
all_IDs2 <- grepl(pattern = "sampled_11_", x = vts_tree$tip.label)
all_IDs2 <- setNames(all_IDs2, vts_tree$tip.label)
toDrop2 <- all_IDs2[all_IDs2 == TRUE]
#to run the infector probability analysis
tree_years <- drop.tip(vts_tree, names(toDrop2))


#get tip names from VirusTreeSimulator tree
tip_names_vts <- tree_years$tip.label

# Change tip names from vts format to ID_migrant format
# change the tips names returned by VirusTreeSimulator to those tip_names
# in the form of ID_migrant
# tip_names1 will have appended to the ID name if it is a migrant or not
# (1 is from region; 2 is from global;
# 12 migrated from region to global;
# 21 migrated from global to region)
# CHECK THIS LINE IN ALL CODES (DELETE THIS LINE LATTER)
tip_names1 <- reorder_tip_names(st_ids_region$tip_name, tip_names_vts)
tip_names1 <- tip_names1[!is.na(tip_names1)]

#to make sure tip names are correct
#compare1 <- unlist(lapply(tip_names1, function(x) str_split(x, "_")[[1]][1]))
#compare2 <- unlist(lapply(tip_names_vts, function(x) str_split(x, "_")[[1]][2]))

#check that 2 vectors are identical
#all(compare1 == compare2)

#get sample number for each sequence
sample_numbers <- unlist(lapply(tip_names_vts,
                                function(x) paste("sample",
                                                  str_split(x, pattern = "_")[[1]][4],
                                                  sep = "_")))



#tree_years$tip.label <- paste(tip_names1, sample_numbers, "count_1", sep = "_")
tree_years$tip.label <- paste("ID", tip_names1, sample_numbers, sep = "_")

# save trees
tree_filename <- paste(prefix_trees, "_sampling.tre", sep="")
write.tree(phy = tree_years, file = tree_filename)


migrant_ID <- unlist(lapply(tree_years$tip.label,
                            function(x) ifelse(str_split(x, "_")[[1]][3] == "1" |
                                               str_split(x, "_")[[1]][3] == "21",
                                               FALSE, TRUE)))
if(sum(migrant_ID) != 0){
  stop("there are mix region and global samples together. Something wrong with
       script to get tip names")
}


#match sampled times to the order of tip names in the phylogenetic tree
sampleTimes <- st_ids_region$sampled_time
sampleTimes <- setNames(sampleTimes, st_ids_region$tip_name)


#Create directory named W (to save everything related to infector probability)
# if it does not exist
if (!dir.exists("output_deepseq/vts/nonRandom")) {
  dir.create("output_deepseq/vts/nonRandom")
}

saveRDS(sampleTimes, paste("output_deepseq/vts/nonRandom/", "sampleTimes.RDS",sep=""))

deep_sequencing_processingnw <- paste("output_deepseq/vts/", "merged_trees_sampling",
                                      "_migrant_years_1_simple_", perc_pop_region,
                                      ".RData", sep="")

save(years, init_sim_date, last_sample_date, start_date,
     end_date_dec, tm, st_ids_region_all, st_ids_region,
     prop_not_region, tree_years, sampleTimes,
     file = deep_sequencing_processingnw)



