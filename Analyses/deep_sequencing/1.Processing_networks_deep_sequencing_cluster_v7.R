# Get transmission matrix and run VirusTreeSimulator
# This code has been tested on a Mac OS and linux and might not work on windows

#Last update was on 28 July 2022

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
simulation_data <- paste(getwd(), "sim_results.tar.gz", sep = "/")
untar(simulation_data)


# This function will generate input file to be used with program
# VirusTreeSimulator and will also run VirusTreeSimulator for each combination of
# inf.csv and sample.csv file.
# It will create a directory "output_deepseq" if directory does not exist
# to save all relevant results.

# Location for VirusTreeSimulator. It should be changed to the correct location
#on your computer.
#Software <- "java -jar VirusTreeSimulator/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
Software <- "java -jar /Applications/VirusTreeSimulator/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
#parameter for VirusTreeSimulator
# parameters for the regional model described in
# Ratman et al. 2016 (Mol Biol Evol 34: 185-203)
# effective population per year growth rate (r) =  2.851904
# -t50 The time point, relative to the time of infection in backwards time, at
# which the population is equal to half its final asymptotic value,
# t50 in Ratman et al paper = -2 years
# -N0 = 0.00593
parameters <- "-demoModel Logistic -N0 0.00593 -growthRate 2.851904 -t50 -2"

#total number of years simulated
years <-  41
#Maximum height to estimate infector probability
#past 5 years
MH <- 5

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

# Read Network data ----
# read departure ID files
# read info from file and convert time from days to decimal years
art_init <- read.csv("results/ART_init.csv")
art_init["time_decimal"] <- days2years(sampleTimes = art_init$time,
                                       init_date = init_sim_date)
dep <- read.csv("results/departure_IDs.csv")
dep["time_decimal"] <- days2years(sampleTimes = dep$time,
                                  init_date = init_sim_date)
names(dep)[2] <- "IDs"

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
#get the average number of people living with HIV at end of simulation
# (last year simulated)
totalPLWHIV <- sum(sim_df$i.num.pop1[tail(sim_df$time, n=365)])/365
PLWHIV <- paste("PLWHIV", totalPLWHIV, sep="_")
#get number of new infections in the past year
newinf_per_year  <- sum(sim_df$incid.pop1[tail(sim_df$time, n=365)])
newinf <- paste("newinf", newinf_per_year, sep="_")

tm <- get_transmat(sim)
tm["year"] <-days2years(tm$at, init_date = init_sim_date)


# Create input files for VirusTreeSimulator ----
#prefix for input files for VirusTreeSimulator
output <- paste("output_deepseq", "output", sep ="/")

#inf file name for VirusTreeSimulator
inf_file <- paste(output, "_inf.csv", sep = "")
#sample file name for VirusTreeSimulator
sample_file <- paste(output, "_sample.csv", sep = "")

#get seed IDs in the network simulations
#seeds are the infected individuals that start a transmission chain
seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

create_inf_csv(tm, time_tr = rep(1980, length(seed_names)), prefix =  output)

# sample IDs and time of sampling
#this sampling strategy will sample IDs of individuals diagnosed, active and not
#on ART
st_ids_region_all <-  sampleIDs(perc = perc_pop_region,
                                start_date = start_date_dec,
                                end_date = end_date_dec,
                                art_init = art_init,
                                departure = dep,
                                diag_info = diag_info,
                                origin = origin,
                                tm = tm,
                                location = "region")


#write proportion of sampled ids that were not in region anymore
prop_not_region <- nrow(st_ids_region_all[(st_ids_region_all$migrant == 2 |
                            st_ids_region_all$migrant == 12), ]) /nrow(st_ids_region_all)

#remove rows in which individuals are not in region anymore
#at time of sampling
st_ids_region <- st_ids_region_all[(st_ids_region_all$migrant == 1 |
                                      st_ids_region_all$migrant == 21), ]

# sampled IDs are not on ART
# Sample 10 viruses per ID to mimic deep-sequencing.
st_ids_region["tip_name"] <-  paste(st_ids_region$sampled_ID,
                                    st_ids_region$migrant,
                                    sep = "_")

#create sample file to be used with VirusTreeSimulator
create_sample_csv2(ids = st_ids_region$sampled_ID,
                   time_seq = st_ids_region$sampled_time,
                   seq_count = 11,
                   prefix = output)


#Create directory named vts (for VirusTreeSimulator) if it does not exist
if (!dir.exists("output_deepseq/vts/")) {
  dir.create("output_deepseq/vts/")
}

#prefix with location of output directory to save results of VirusTreeSimulator
prefix_vts <-  "output_deepseq/vts/results_vts"

# Run VirusTreeSimulator ----
cmd <-  paste(Software, parameters, inf_file, sample_file, prefix_vts, sep = " ")
system(cmd)


#Create directory to save phylogenetic trees
if (!dir.exists("output_deepseq/vts/merged_trees/")) {
  dir.create("output_deepseq/vts/merged_trees")
}

#prefix with location of output directory to save merged phylogenetic trees
prefix_trees <-  "output_deepseq/vts/merged_trees/results"

#Analyse trees generated with VirusTreeSimulator ----
#read all trees from vts
list_trees <- dir("output_deepseq/vts",
                  pattern = "*_detailed.nex",
                  full.names = TRUE)
trees <- lapply(list_trees, read.nexus)

#add root.edge
if(length(trees) > 1){
  trees_rootedge <- lapply(trees, add_root_edge2,
                           root.edge_value = 0)
  vts_tree <- merge_trees(trees_rootedge)
} else{
  trees_rootedge <- add_root_edge2(tree = trees[[1]], root.edge_value = 0)
  vts_tree <- trees_rootedge
}

#now keep only one sequence per ID
#to create a consensus tree and
#to filter tree by infector probability
all_IDs <- grepl(pattern = "sampled_10", x = vts_tree$tip.label)
all_IDs <- setNames(all_IDs, vts_tree$tip.label)
toDrop <- all_IDs[all_IDs == FALSE]

#create a consesus sequences by keeping just one virus per ID
consensus_tree <- drop.tip(vts_tree, names(toDrop))

# Change tip names from vts format to ID_migrant format
# change the tips names returned by VirusTreeSimulator to those tip_names
# in the form of ID_migrant
tip_names_migrant <- st_ids_region$tip_name
consensus_tree$tip.label <- reorder_tip_names(tip_names_migrant,
                                              consensus_tree$tip.label)


#now remove only sampled_11. I simulated 11 sequences per ID
#because the detailed.nex output of the VirusTreeSimulator
#when read into R generated a weird number of internal numbers
#most likely because of unifurcations (read it here https://github.com/niemasd/CoaTran)
#By removing some sequences from the tree, the number of internal nodes are then
#fixed to the correct number.
all_IDs2 <- grepl(pattern = "sampled_11_", x = vts_tree$tip.label)
all_IDs2 <- setNames(all_IDs2, vts_tree$tip.label)
toDrop2 <- all_IDs2[all_IDs2 == TRUE]

#this the final tree in which branch lengths are in unit of calendar time
tree_years <- drop.tip(vts_tree, names(toDrop2))


#get tip names from VirusTreeSimulator tree
tip_names_vts <- tree_years$tip.label

# Change tip names from vts format to ID_migrant format
# change the tips names returned by VirusTreeSimulator to those st_ids_region$tip_name
# in the form of ID_migrant
# tip_names1 will have appended to the ID name if it is a migrant or not
# (1 is from region; 2 is from global;
# 12 migrated from region to global;
# 21 migrated from global to region)
tip_names1 <- reorder_tip_names(st_ids_region$tip_name, tip_names_vts)
tip_names1 <- tip_names1[!is.na(tip_names1)]


#get sample number for each sequence
sample_numbers <- unlist(lapply(tip_names_vts,
                                function(x) paste("sample",
                                                  str_split(x, pattern = "_")[[1]][4],
                                                  sep = "_")))


tree_years$tip.label <-  paste("ID", tip_names1, sample_numbers, sep = "_")

# save trees
tree_filename <- paste(prefix_trees, "_sampling.tre", sep = "")
write.tree(phy = tree_years, file = tree_filename)

#save consensus sequence
tree_filename_consensus <- paste(prefix_trees, "_sampling_consensus.tre", sep = "")
write.tree(phy = consensus_tree, file = tree_filename_consensus)

migrant_ID <- unlist(lapply(tree_years$tip.label,
                            function(x)
                              ifelse(
                                str_split(x, "_")[[1]][3] == "1" |
                                  str_split(x, "_")[[1]][3] == "21",
                                FALSE,
                                TRUE
                              )))
if (sum(migrant_ID) != 0) {
  stop(
    "there are mix region and global samples together. Something wrong with
         script to get tip names"
  )
}

#Prepare data to calculate infector probabilities ----
all_cd4s <- get_cd4s_sampling(st_ids_region, stages)

tips <- unlist(lapply(consensus_tree$tip.label, function(x) str_split(x, "_")[[1]][1]))

#match cd4s to order of tip names in the phylogenetic tree
all_cd4s <- all_cd4s[match(tips, names(all_cd4s))]
all_cd4s <- setNames(all_cd4s, consensus_tree$tip.label)

#get ehi (early HIV infection)
#get which ids have CD4 count equivalent to stage 0 (Early HIV stage
#of infection) for the past 6 months at time of sampling
recent_ids <- recency_test(st_ids_region, stages)

#named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
ehis <- ifelse(names(all_cd4s) %in% recent_ids, TRUE, FALSE)
ehis <- setNames(ehis, names(all_cd4s))

#match sampled times to the order of tip names in the phylogenetic tree
sampleTimes <- st_ids_region$sampled_time
sampleTimes <- setNames(sampleTimes, st_ids_region$tip_name)


#Calculate infector probabilities ----

# to calculate infector probabilities
# sampleTimes: must use years
# cd4s: named numeric vector, cd4 at time of sampling
# ehi: named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
# numberPeopleLivingWithHIV: scalar
# numberNewInfectionsPerYear: scalar

#only infector probabilities with MH will be calculated
W <- phylo.source.attribution.hiv.msm(consensus_tree,
                                      sampleTimes[consensus_tree$tip.label],
                                      cd4s = all_cd4s[consensus_tree$tip.label],
                                      ehi = ehis[consensus_tree$tip.label],
                                      numberPeopleLivingWithHIV  = totalPLWHIV,
                                      numberNewInfectionsPerYear = newinf_per_year,
                                      maxHeight = MH,
                                      res = 1e3,
                                      treeErrorTol = Inf)


#convert list to dataframe
W1 <- as.data.frame(W)



W1["donor_ID"] <- unlist(lapply(W$donor, function(x) str_split(x, pattern = "_")[[1]][1]))
W1["recip_ID"] <- unlist(lapply(W$recip, function(x) str_split(x, pattern = "_")[[1]][1]))
W1$donor_ID <- as.numeric(W1$donor_ID)
W1$recip_ID <- as.numeric(W1$recip_ID)

#subset data
Wsub <- data.frame(donor_ID = W1$donor_ID,
                   recip_ID = W1$recip_ID,
                   infectorProbability = W1$infectorProbability)

Wsub$donor_ID <- as.integer(Wsub$donor_ID)
Wsub$recip_ID <- as.integer(Wsub$recip_ID)


#get individuals that seroconverted in region within the transmission matrix
tm_region <-  subset(tm, infOrigin == "region" & susOrigin == "region")

#keep_only_rows that tips are in phylogenetic tree
rows_to_keep <- keep_row(df = tm_region, tree = consensus_tree)

if (!is.null(rows_to_keep)) {
  tm1 <- tm_region[rows_to_keep, ]
} else{
  #if there is no rows to keep, assign all values of tm1 to NA
  tm1 <- tm_region[1, ]
  tm1[1, ] <- NA
}

#create a dataframe similar to the Wsub to use with function semi_join
tm_all1 <- data.frame(donor_ID = tm1$inf,
                      recip_ID = tm1$sus,
                      infectorProbability = 1)

#get real transmissions from donor (inf) to recipient (sus)
real_trans <- semi_join(Wsub, tm_all1, by = c("donor_ID", "recip_ID"))


#Create directory named W (to save everything related to infector probability)
# if it does not exist
if (!dir.exists("output_deepseq/vts/W")) {
  dir.create("output_deepseq/vts/W")
}

saveRDS(sampleTimes,
        paste("output_deepseq/vts/W/", "sampleTimes.RDS", sep = ""))

deep_sequencing_processingnw <- paste("output_deepseq/vts/",
                                      "merged_trees_sampling",
                                      "_migrant_years_1_simple_",
                                      perc_pop_region,
                                      ".RData",
                                      sep = "")

save(years,
     MH,
     init_sim_date,
     last_sample_date,
     start_date,
     end_date_dec,
     tm,
     st_ids_region,
     prop_not_region,
     consensus_tree,
     tree_years,
     sampleTimes,
     all_cd4s,
     ehis,
     newinf_per_year,
     totalPLWHIV,
     W,
     W1,
     tm_all1,
     real_trans,
     file = deep_sequencing_processingnw)

summaryW(sim = perc_pop_region,
         tm = tm,
         W,
         tree = consensus_tree,
         code = "TrueTrees",
         prefix = NULL,
         labels = TRUE)



#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time

processing_network_time <- data.frame(start = start_time, end = end_time)
saveRDS(processing_network_time, "processing_network_time_deepseq.RDS")

