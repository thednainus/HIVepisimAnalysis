# 05 May 2022
# Get transmission matrix ----
# and run VirusTreeSimulator
# This code has been tested on a Mac OS and might not work on windows or linux

# I also only analyse the firt 30 replicates
# In this script I use a vector for incidence and prevalece for calculation of
#infector probabilities

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
line_number <- as.numeric(commandArgs(trailingOnly = TRUE))

params <- readRDS(system.file("data/simulations_imperial_cluster.RDS",
                              package = "HIVepisimAnalysis"))

params <- params[c(1:30, 101:130, 201:230, 301:330, 401:430,
                   501:530, 601:630, 701:730, 801:830, 901:930,
                   1001:1030, 1101:1130, 1201:1230, 1301:1330, 1401:1430,
                   1501:1530, 1601:1630, 1701:1730, 1801:1830, 1901:1930),]

#params <- params[c(3:30, 101:130, 201:230, 301:330, 401:430,
#                   501:530, 601:630, 701:730, 801:830, 901:930,
#                   1001:1030, 1101:1130, 1201:1230, 1301:1330, 1401:1430,
#                   1501:1530, 1601:1630, 1701:1730, 1801:1830, 1901:1930),]

params <- params[line_number,]


perc_pop_region <- params$perc
#perc_pop_global <- 3 * perc_pop_region
perc_pop_global <- 0.05


#untar tar file
simulation_dir <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants/params_1067/rep_28"
#simulation_dir <- paste("/rds/general/user/fferre15/ephemeral/",
#                         paste("params", params$params, sep = "_"),
#                         "/",
#                         paste("rep", params$rep, sep = "_"), sep = "")
simulation_data <- paste(simulation_dir, "sim_results.tar.gz", sep = "/")

mkdir_info <- paste(simulation_dir,
                    paste("perc", params$perc, sep = "_"),
                    sep = "/")
where2save <- paste(mkdir_info, ".", sep = "/")

write.table(mkdir_info, file = "mkdir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)

write.table(where2save, file = "dir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)

untar(simulation_data)




# This function will generate input file to be used with program
# VirusTreeSimulator
# and will also run VirusTreeSimulator for each combination of
# inf and sample file
# It will create a directory "output" if directory does not exist
# and will save results to the directory "output"

# Location for VirusTreeSimulator. It should be changed to the correct location on your computer.
#Software <- "java -jar /Applications/VirusTreeSimulator/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#Software <- "java -jar VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
#parameters <- "-demoModel Constant -N0 1"
# parameters following Ratman et al. 2016 (Mol Biol Evol 34: 185-203)
# effective population per year growth rate (r) =  2.851904
#  -t50 The time point, relative to the time of infection in backwards time, at
# which the population is equal to half its final asymptotic value,
# t50 in Ratman et al paper = -2 years
parameters <- "-demoModel Logistic -N0 1 -growthRate 2.851904 -t50 -2"

#maximum height
#analyse the past 15 years
#MH <-  15
# total number of years simulated
years <-  41
#area and max_value is used in function get_tipNames
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
#sampling within the last 15 years
start_date <- ymd("2006-01-01")
start_date_dec <- decimal_date(start_date)
end_date_dec <- decimal_date(last_sample_date)



#Create directory named output if it does not exist
if (!dir.exists("output")) {
  dir.create("output")
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


sim <- readRDS("results/results_sim.RDS")
sim_df <- as.data.frame(sim)
#convert time to years decimal
sim_df["years"] <- days2years(sampleTimes = sim_df$time, init_date = init_sim_date)
#convert days to years for tree branch lengths
sim_df["years_branch_lengths"] <- sim_df$time * 1/365


#People living with HIV (PLWHIV) in the past year
totalPLWHIV <- sum(sim_df$i.num.pop1[tail(sim_df$time, n=365)])/365
PLWHIV <- paste("PLWHIV", totalPLWHIV, sep="_")



#get number of new infections in the past year
newinf_per_year  <- sum(sim_df$incid.pop1[tail(sim_df$time, n=365)])
newinf <- paste("newinf", newinf_per_year, sep="_")

#get transmat
tm <- get_transmat(sim)

# get transmat by seed
# here I only get the ones that started in region
# later check all trees
#tm2 <- get.transmat.phylo(tm, by_areas = area, max_value = max_value)
#seed_names <- names(tm2)
# if using the list of transmission matrix by seed
if(!is.null(tm)){

    # Get tip names in the form of ID_migrant
    tip_names <- get_tipNames(tm,
                              format = "migrant",
                              by_areas = area,
                              max_value = max_value)

    # check number of individual in network to make sure it is not empty
    tip_names_migrant_ID <- unlist(lapply(tip_names,
                                          function(x) str_split(x, "_")[[1]][2]))


    if(length(tip_names_migrant_ID) > 0){


      output <- paste("output", "output", sep ="/")
      output <- paste(output, PLWHIV, newinf, sep = "_")

      #inf file name for VirusTreeSimulator
      inf_file <- paste(output, "_inf.csv", sep = "")
      #sample file name for VirusTreeSimulator
      sample_file <- paste(output, "_sample.csv", sep = "")

      seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

      create_inf_csv(tm, time_tr = rep(0, length(seed_names)), prefix = output)

      # sample IDs from region and time of sampling
      st_ids_region_all <- sampleIDs2(perc = perc_pop_region, start_date = start_date_dec,
                                      end_date = end_date_dec, art_init = art_init,
                                      departure = dep, diag_info = diag_info,
                                      origin = origin,
                                      tm = tm, location = "region")



      #write proportion of sampled ids that were not in region anymore
      prop_not_region <- nrow(st_ids_region_all[(st_ids_region_all$migrant == 2 |
                                                st_ids_region_all$migrant == 12),])/nrow(st_ids_region_all)

      #remove rows in which individuals are not in region anymore
      #before or at time of sampling
      st_ids_region <- st_ids_region_all[(st_ids_region_all$migrant == 1 |
                                          st_ids_region_all$migrant == 21 ),]

      # sampled IDs are not on ART

      st_ids_region["date"] <- as.Date(date_decimal(st_ids_region$sampled_time))
      st_ids_region["time_days"] <- as.numeric(st_ids_region$date - init_sim_date)
      st_ids_region["time_years"] <- st_ids_region$time_days * (1/365)
      st_ids_region["tip_name"] <- paste(st_ids_region$sampled_ID, st_ids_region$migrant,
                                         sep = "_")


      # sample IDs from global and time of sampling
      #here does not really matter if some individuals have migrated to region
      #because this data will be used just to set the clock rate
      #so individuals infected at the begining of the epidemic can be sampled
      st_ids_global_all <- sampleIDs2(perc = perc_pop_global, start_date = decimal_date(init_sim_date),
                                 end_date = decimal_date(last_sample_date), art_init = art_init,
                                 departure = dep, diag_info = diag_info,
                                 origin = origin,
                                 tm = tm, location = "global")


      #remove rows in which individuals are not in global anymore
      #before or at time of sampling
      st_ids_global <- st_ids_global_all[(st_ids_global_all$migrant == 2 |
                                          st_ids_global_all$migrant == 12 ),]


      st_ids_global["date"] <- as.Date(date_decimal(st_ids_global$sampled_time))
      st_ids_global["time_days"] <- as.numeric(st_ids_global$date - init_sim_date)
      st_ids_global["time_years"] <- st_ids_global$time_days * (1/365)
      st_ids_global["tip_name"] <- paste(st_ids_global$sampled_ID, st_ids_global$migrant,
                                         sep = "_")

      #get stage of HIV infection at time of sampling for each individual


      # here I had to simulate 2 sequences per ID because there is no tree with
      # 1 individual.
      # Because I am applying a sampling process, I may have trees with 1 individual
      # and the following scripts would not work
      create_sample_csv2(ids = c(st_ids_region$sampled_ID, st_ids_global$sampled_ID),
                         time_seq = c(st_ids_region$time_years, st_ids_global$time_years),
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
  list_trees <- dir("output/vts", pattern = "*_detailed.nex", full.names = TRUE)
  trees <- lapply(list_trees, read.nexus)
  #add root.edge

  #trees_rootedge <- lapply(trees, add_root_edge,
  #                         total_sim_steps = decimal_date(end_sim_date) - decimal_date(init_sim_date),
  #                         root.edge_value = 0)
  trees_rootedge <- lapply(trees, add_root_edge2,
                           root.edge_value = 0)
  vts_tree <- merge_trees(trees_rootedge)

  #now keep only one sequence per ID
  #drop tips in which node has died before end of simulation
  all_IDs <- grepl(pattern = "sampled_2", x = vts_tree$tip.label)
  all_IDs <- setNames(all_IDs, vts_tree$tip.label)
  toDrop <- all_IDs[all_IDs == TRUE]
  consensus_tree <- drop.tip(vts_tree, names(toDrop))


  # convert branch lengths from days to years
  #tree_years <- convert_branches(tree = consensus_tree, scale = 1/365)
  tree_years <- consensus_tree

  #get tip names from VirusTreeSimulator tree
  tip_names_vts <- tree_years$tip.label

  # Change tip names from vts format to ID_migrant format
  # change the tips names returned by VirusTreeSimulator to those tip_names
  # in the form of ID_migrant
  tip_names_migrant <- c(st_ids_region$tip_name, st_ids_global$tip_name)
  tree_years$tip.label <- reorder_tip_names(tip_names_migrant, tip_names_vts)


  # save tree to simulate sequence alignment using Python script

  tree_filename <- paste(prefix_vts, "_merged_trees_sampling_",
                         perc_pop_region, ".tre", sep="")
  write.tree(phy = tree_years, file = tree_filename)


  # Calculate infector probability ----
  all_cd4s <- get_cd4s_sampling(rbind(st_ids_region, st_ids_global), stages)


  tips <- unlist(lapply(tree_years$tip.label, function(x) str_split(x, "_")[[1]][1]))

  #match cd4s to order of tip names in the phylogenetic tree
  all_cd4s <- all_cd4s[match(tips, names(all_cd4s))]
  all_cd4s <- setNames(all_cd4s, tree_years$tip.label)

  #get ehi (early HIV infection)
  #get which ids that have CD4 count equivalent to stage 0 (Early HIV stage
  #of infection) is at this stage for the past 6 months at time of sampling
  recent_ids <- recency_test(rbind(st_ids_region, st_ids_global), stages)
  #which(names(all_cd4s) %in% recent_ids)

  #ehis = early hiv stage of infection
  #named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )

  #ehis <- ifelse(all_cd4s == 1e3, TRUE, FALSE)
  ehis <- ifelse(names(all_cd4s) %in% recent_ids, TRUE, FALSE)
  ehis <- setNames(ehis, names(all_cd4s))

  #match sampled times to the order of tip names in the phylogenetic tree
  sampleTimes <- c(st_ids_region$sampled_time, st_ids_global$sampled_time)
  sampleTimes <- setNames(sampleTimes, c(st_ids_region$tip_name, st_ids_global$tip_name))


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



  W <- phylo.source.attribution.hiv.msm(onlyregion_tree, sampleTimes[onlyregion_tree$tip.label],
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

  W_filename <- paste("output/vts/W/", "merged_trees_sampling", "_migrant_years_1_simple_",
                      perc_pop_region, ".RData", sep="")
  save(years, max_value, init_sim_date, last_sample_date, start_date,
       end_date_dec, tm, st_ids_region_all, st_ids_region, st_ids_global_all,
       st_ids_global, prop_not_region, tree_years, onlyregion_tree,
       sampleTimes, all_cd4s, ehis, newinf_per_year, totalPLWHIV, W,
       file = W_filename)


  summaryW(sim = perc_pop_region, tm = tm, W,
           tree = onlyregion_tree, code = "TrueTrees",
           prefix = NULL, labels = TRUE)

}


#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time

processing_network_time <- data.frame(start = start_time, end = end_time)
saveRDS(processing_network_time, "processing_network_time.RDS")

