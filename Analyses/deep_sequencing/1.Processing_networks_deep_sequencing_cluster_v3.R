# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux

# this script will add a function to sample IDs and sample times

#start of script
start_time <- Sys.time()

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


# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output" if directory does not exist
# and will save results to the directory "output"

# Location for VirusTreeSimulator. It should be changed to the correct location on your computer.
#Software <- "java -jar VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
#parameters <- "-demoModel Constant -N0 1"
# parameters following Ratman et al. 2016 (Mol Biol Evol 34: 185-203)
# parameters in Ratman et al. 2016 is per year, and below I converted it to
# days as my simulations are in units of days
# effective population per year growth rate (r) =  2.851904
# then effective population per day growth rate (r) =  0.007813436
#  -t50 The time point, relative to the time of infection in backwards time, at
# which the population is equal to half its final asymptotic value,
# t50 in Ratman et al paper = -2 years
# t50 in days = 2/(2*365) = 0.002739726
parameters <- "-demoModel Logistic -N0 1 -growthRate 0.007813436 -t50 -0.002739726"

years <-  40
area <-  "all"
max_value <-  NULL
#beginning of simulation time
init_sim_date <- ymd("1980-01-01")
# year of last sample dates in simulations
last_sample_date <- as.Date(x = years*365, origin = init_sim_date)


#Times for sampling IDs and sampling times
start_date <- ymd("1996-01-01")
start_date_dec <- decimal_date(start_date)
end_date <- ymd("2015-06-30")
end_date_dec <- decimal_date(end_date)


#Create directory named output_deepseq if it does not exist
if (!dir.exists("output_deepseq")) {
  dir.create("output_deepseq")
}

# read departure ID files
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

#read simulation results
sim <- readRDS("run1/results_sim.RDS")
sim_df <- as.data.frame(sim)

tm <- get_transmat(sim)


if(!is.null(tm)){

    # Get tip names in the form of ID_migrant
    tip_names <- get_tipNames(tm, format = "migrant",
                                   by_areas = area, max_value = max_value)

    # check number of individual within "region"
    # region is code as 1 and 21
    tip_names_migrant_ID <- unlist(lapply(tip_names, function(x) str_split(x, "_")[[1]][2]))
    #total_tips_region <- sum(tip_names_migrant_ID == "1" | tip_names_migrant_ID == "21")

    if(length(tip_names_migrant_ID) > 0){


      output <- paste("output_deepseq", "output", sep ="/")
      #output <- paste(output, PLWHIV, newinf, sep = "_")

      #inf file name for VirusTreeSimulator
      inf_file <- paste(output, "_inf.csv", sep = "")
      #sample file name for VirusTreeSimulator
      sample_file <- paste(output, "_sample.csv", sep = "")

      seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

      create_inf_csv(tm, time_tr = rep(0, length(seed_names)), prefix=output)

      # sample IDs and time of sampling
      st_ids_region <- sampleIDs(perc = 0.05, start_date = start_date_dec,
                                 end_date = end_date_dec, art_init = art_init,
                                 departure = dep, diag_info = diag_info,
                                 tm = tm, location = "region")

      # sampled IDs are not on ART
      # Sample 10 viruses per ID
      # note that this is just to start the simulations
      # but I should try to find real values for those

      st_ids_region["date"] <- as.Date(date_decimal(st_ids_region$sampled_time))
      st_ids_region["time_days"] <- as.numeric(st_ids_region$date - init_sim_date)

      create_sample_csv2(ids = st_ids_region$sampled_ID,
                         time_seq = st_ids_region$time_days,
                         seq_count = 10, prefix = output)

      #Create directory named VTS (for VirusTreeSimulator) if it does not exist
      if (!dir.exists("output_deepseq/vts/")) {
        dir.create("output_deepseq/vts/")
      }

      #prefix with location of output directory to save results of VirusTreeSimulator
      prefix_vts <-  "output_deepseq/vts/results_vts"

      # Run VirusTreeSimulator
      cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
      system(cmd)
    }

  #Create directory to save phylogenetic trees
  if (!dir.exists("output_deepseq/vts/merged_trees/")) {
    dir.create("output_deepseq/vts/merged_trees")
  }

  #prefix with location of output directory to save merged phylogenetic trees
  prefix_trees <-  "output_deepseq/vts/merged_trees/results"

  #read all trees from vts
  list_trees <- dir("output_deepseq/vts", pattern = "*_simple.nex", full.names = TRUE)
  trees <- lapply(list_trees, read.nexus)
  #add root.edge
  trees_rootedge <- lapply(trees, add_root_edge, total_sim_steps = years * 365,
                           root.edge_value = 0)
  vts_tree <- merge_trees(trees_rootedge)


  # convert branch lengths from days to years
  tree_years <- convert_branches(tree = vts_tree, scale = 1/365)


  #get tip names from VirusTreeSimulator tree
  tip_names_vts <- tree_years$tip.label

  # Change tip names from vts format to ID_migrant format
  # change the tips names returned by VirusTreeSimulator to those tip_names
  # in the form of ID_migrant
  # tip_names1 will have appended to the ID name if it is a migrant or not
  # (1 is from region; 2 is from global;
  # 12 migrated from region to global;
  # 21 migrated from global to region)
  tip_names1 <- reorder_tip_names(tip_names, tip_names_vts)
  tip_names1 <- tip_names1[!is.na(tip_names1)]

  #get sample number for each sequence
  sample_numbers <- unlist(lapply(tip_names_vts, function(x) paste("sample",
                                                                   str_split(x, pattern = "_")[[1]][4],
                                                                   sep = "_")
    ))


  #tree_years$tip.label <- paste(tip_names1, sample_numbers, "count_1", sep = "_")
  tree_years$tip.label <- paste("ID", tip_names1, sample_numbers, sep = "_")

  # save trees
  tree_filename <- paste(prefix_trees, "_sampling.tre", sep="")
  write.tree(phy = tree_years, file = tree_filename)

  # Remove sequences from "global" -----
  # to calculate infector probability
  migrant_ID <- unlist(lapply(tree_years$tip.label,
                              function(x) ifelse(str_split(x, "_")[[1]][3] == "1" | str_split(x, "_")[[1]][3] == "21",
                                                 FALSE, TRUE)))
  if(sum(migrant_ID) != 0){
    stop("there are mix region and global samples together. Something wrong with
         script to get tip names")
  }



}

#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time

processing_network_time <- data.frame(start = start_time, end = end_time)
saveRDS(processing_network_time, "processing_network_time.RDS")

