# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux

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

# percentage of population to sampled IDs
#perc_pop_region <- commandArgs(trailingOnly = TRUE)
perc_pop_region <- 0.05
perc_pop_global <- 3 * perc_pop_region






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
# year of last sample dates in simulations
#beginning of simulation time
#first case of HIV in San Diego was reported in 1980
init_sim_date <- ymd("1980-01-01")
# end of simulation date
end_sim_date <- ymd("2021-01-01")
# year of last sample dates in simulations
last_sample_date <- end_sim_date


#Times for sampling IDs and sampling times for within region
start_date <- ymd("1996-01-01")
start_date_dec <- decimal_date(start_date)
end_date <- ymd("2015-06-30")
end_date_dec <- decimal_date(end_date)



#Create directory named output if it does not exist
if (!dir.exists("output")) {
  dir.create("output")
}


# read departure ID files
#read info from file and convert time from days to decimal years
art_init <- read.csv("ART_init.csv")
art_init["time_decimal"] <- days2years(sampleTimes = art_init$time,
                                       init_date = init_sim_date)
dep <- read.csv("departure_IDs.csv")
dep["time_decimal"] <- days2years(sampleTimes = dep$time,
                                  init_date = init_sim_date)

diag_info <- read.csv("diag_time.csv")
diag_info["time_decimal"] <- days2years(sampleTimes = diag_info$time,
                                        init_date = init_sim_date)
stages <- read.csv("stages.csv")
stages["time_decimal"] <- days2years(sampleTimes = stages$time,
                                     init_date = init_sim_date)


sim <- readRDS("results_sim.RDS")
sim_df <- as.data.frame(sim)
#convert time to years decimal
sim_df["years"] <- days2years(sampleTimes = sim_df$time, init_date = init_sim_date)
#get the average number of people living with HIV at end of simulation
# (last year simulated)
#get total people living with HIV relative to the sampling times for region
# end_date_dec = 2015.493
total_days <- length(sim_df$i.num.pop1[sim_df$years >= 2015 & sim_df$years < 2016])
totalPLWHIV <- sum(sim_df$i.num.pop1[sim_df$years >= 2015 & sim_df$years < 2016])/total_days
#totalPLWHIV <- sum(sim_df$i.num.pop1[tail(sim_df$time, n=365)])/365
PLWHIV <- paste("PLWHIV", totalPLWHIV, sep="_")
#get number of new infections in the past year
#get number of new infections relative to the sampling times for region
# end_date_dec = 2015.493
newinf_per_year <- sum(sim_df$incid.pop1[sim_df$years >= 2015 & sim_df$years < 2016])/total_days
#newinf_per_year  <- sum(sim_df$incid.pop1[tail(sim_df$time, n=365)])
newinf <- paste("newinf", newinf_per_year, sep="_")

tm <- get_transmat(sim)

# get transmat by seed
# here I only get the ones that started in region
# later check all trees
#tm2 <- get.transmat.phylo(tm, by_areas = area, max_value = max_value)
#seed_names <- names(tm2)
# if using the list of transmission matrix by seed
if(!is.null(tm)){

    # Get tip names in the form of ID_migrant
    tip_names <- get_tipNames(tm, format = "migrant",
                                   by_areas = area, max_value = max_value)

    # check number of individual within "region"
    # region is code as 1 and 21
    tip_names_migrant_ID <- unlist(lapply(tip_names, function(x) str_split(x, "_")[[1]][2]))
    #total_tips_region <- sum(tip_names_migrant_ID == "1" | tip_names_migrant_ID == "21")

    if(length(tip_names_migrant_ID) > 0){


      output <- paste("output", "output", sep ="/")
      output <- paste(output, PLWHIV, newinf, sep = "_")

      #inf file name for VirusTreeSimulator
      inf_file <- paste(output, "_inf.csv", sep = "")
      #sample file name for VirusTreeSimulator
      sample_file <- paste(output, "_sample.csv", sep = "")

      seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

      create_inf_csv(tm, time_tr = rep(0, length(seed_names)), prefix=output)

      # sample IDs from region and time of sampling
      st_ids_region <- sampleIDs(perc = perc_pop_region, start_date = start_date_dec,
                                 end_date = end_date_dec, art_init = art_init,
                                 departure = dep, diag_info = diag_info,
                                 tm = tm, location = "region")

      # sampled IDs are not on ART

      st_ids_region["date"] <- as.Date(date_decimal(st_ids_region$sampled_time))
      st_ids_region["time_days"] <- as.numeric(st_ids_region$date - init_sim_date)


      # sample IDs from global and time of sampling
      st_ids_global <- sampleIDs(perc = perc_pop_global, start_date = decimal_date(init_sim_date),
                                 end_date = decimal_date(last_sample_date), art_init = art_init,
                                 departure = dep, diag_info = diag_info,
                                 tm = tm, location = "global")

      # sampled IDs are not on ART

      st_ids_global["date"] <- as.Date(date_decimal(st_ids_global$sampled_time))
      st_ids_global["time_days"] <- as.numeric(st_ids_global$date - init_sim_date)

      #get stage of HIV infection at time of sampling for each individual


      # here I had to simulate 2 sequences per ID because there is not tree with
      # 1 individual.
      # Because I am applying a sampling process, I may have trees with 1 individual
      # and the following scripts would not work
      create_sample_csv2(ids = c(st_ids_region$sampled_ID, st_ids_global$sampled_ID),
                         time_seq = c(st_ids_region$time_days, st_ids_global$time_days),
                         seq_count = 2, prefix = output)


      #Create directory named VTS (for VirusTreeSimulator) if it does not exist
      if (!dir.exists("output/vts/")) {
        dir.create("output/vts/")
      }

      #prefix with location of output directory to save results of VirusTreeSimulator
      prefix_vts <-  "output/vts/results_vts"

      # Run VirusTreeSimulator
      cmd <- paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
      system(cmd)
    }

  #read all trees from vts
  list_trees <- dir("output/vts", pattern = "*_simple.nex", full.names = TRUE)
  trees <- lapply(list_trees, read.nexus)
  #add root.edge
  trees_rootedge <- lapply(trees, add_root_edge, total_sim_steps = 365 * years)
  vts_tree <- merge_trees(trees_rootedge)

  #now keep only one sequence per ID
  #drop tips in which node has died before end of simulation
  all_IDs <- grepl(pattern = "sampled_2", x = vts_tree$tip.label)
  all_IDs <- setNames(all_IDs, vts_tree$tip.label)
  toDrop <- all_IDs[all_IDs == TRUE]
  consensus_tree <- drop.tip(vts_tree, names(toDrop))


  # convert branch lengths from days to years
  tree_years <- convert_branches(tree = consensus_tree, scale = 1/365)

  #get tip names from VirusTreeSimulato tree
  tip_names_vts <- tree_years$tip.label

  # Change tip names from vts format to ID_migrant format
  # change the tips names returned by VirusTreeSimulator to those tip_names
  # in the form of ID_migrant
  tree_years$tip.label <- reorder_tip_names(tip_names, tip_names_vts)


  # save tree to simulate sequence alignment using Python script

  tree_filename <- paste(prefix_vts, "_merged_trees_sampling_", perc_pop_region, ".tre", sep="")
  write.tree(phy = tree_years, file = tree_filename)


  # Calculate infector probability ----
  all_cd4s <- get_cd4s_sampling(rbind(st_ids_region, st_ids_global), stages)

  tips <- unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][1]))

  #match cd4s to order of tip names in the phylogenetic tree
  all_cd4s <- all_cd4s[match(tips, names(all_cd4s))]
  all_cd4s <- setNames(all_cd4s, tree_years$tip.label)

  #get ehi (early HIV infection)
  #named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  ehis <- ifelse(all_cd4s == 1e3, TRUE, FALSE)

  #match sampled times to the order of tip names in the phylogenetic tree
  sampleTimes <- c(st_ids_region$sampled_time, st_ids_global$sampled_time)
  sampleTimes <- setNames(sampleTimes, c(st_ids_region$sampled_ID, st_ids_global$sampled_ID))
  sampleTimes <- sampleTimes[match(tips, names(sampleTimes))]
  sampleTimes <- setNames(sampleTimes, tree_years$tip.label)


  # to calculate infector probabilities
  # sampleTimes: must use years
  # cd4s: named numeric vector, cd4 at time of sampling
  # ehi: named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  # numberPeopleLivingWithHIV: scalar
  # numberNewInfectionsPerYear: scalar

  #drop tips from global to calculate infector probability

  tipLocation <- grepl(pattern = "_2$|_12", x = tree_years$tip.label, perl = TRUE)
  tipLocation <- setNames(tipLocation, tree_years$tip.label)
  tips2drop <- tipLocation[tipLocation == TRUE]

  onlyregion_tree <- drop.tip(tree_years, names(tips2drop))

  W <- phylo.source.attribution.hiv.msm( onlyregion_tree, sampleTimes[onlyregion_tree$tip.label],
                                         cd4s = all_cd4s[onlyregion_tree$tip.label],
                                         ehi = ehis[onlyregion_tree$tip.label],
                                         numberPeopleLivingWithHIV  = totalPLWHIV,
                                         numberNewInfectionsPerYear = newinf_per_year,
                                         maxHeight = years,
                                         res = 1e3,
                                         treeErrorTol = Inf)

  #Create directory named W (to save everything related to infector probability)
  # if it does not exist
  if (!dir.exists("output/vts/W")) {
    dir.create("output/vts/W")
  }

  saveRDS(sampleTimes, paste("output/vts/W/", "sampleTimes.RDS",sep=""))
  #prefix <- "test"

  W_filename <- paste("output/vts/W/", "merged_trees_sampling", "_migrant_years_1_simple_", perc_pop_region, ".RData", sep="")
  save(years, MH = years, max_value, init_sim_date, last_sample_date, start_date,
       end_date, tm, st_ids_region, st_ids_global,
       tree_years, sampleTimes, all_cd4s, ehis, newinf_per_year, totalPLWHIV, W,
       file = W_filename)

  summaryW(sim = "1", tm = tm, W, ID = "sampled",
            tree = tree_years, code = "TrueTrees",
           prefix = NULL, labels = TRUE)

}


#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time

processing_network_time <- data.frame(start = start_time, end = end_time)
saveRDS(processing_network_time, "processing_network_time.RDS")
