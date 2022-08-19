# Estimate phylogenetic tree using IQ-TREE
# this script takes as argument "1000" or "10000"
# to specify alignments of 1000bp or 10000bp
# start of script

#In this script I am analysing 30 replicates per combination of parameter values.
start_time <- Sys.time()


# and convert branch lengths to unit of calendar time
library(treedater)
library(DescTools)
library(phydynR)
library(ape)
library(stringr)
library(HIVepisimAnalysis)
library(dplyr)

# percentage of population to sampled IDs
get_params <- commandArgs(trailingOnly = TRUE)

line_number <-as.numeric(get_params[1])

params <- readRDS(system.file("data/simulations_imperial_cluster.RDS",
                              package = "HIVepisimAnalysis"))

params <- params[c(1:30, 101:130, 201:230, 301:330, 401:430,
                   501:530, 601:630, 701:730, 801:830, 901:930,
                   1001:1030, 1101:1130, 1201:1230, 1301:1330, 1401:1430,
                   1501:1530, 1601:1630, 1701:1730, 1801:1830, 1901:1930),]

params <- params[line_number,]

seqlength <- as.numeric(get_params[2])
seq_length <- paste(seqlength, "bp", sep = "")

#get where the data is located
migration_variation <- as.character(get_params[3])
#migration_variation <- paste("best_trajectories", migration_variation, sep ="_")
simulation_dir <- paste("/rds/general/user/fferre15/ephemeral/",
                        migration_variation, "/",
                        paste("params", params$params, sep = "_"),
                        "/",
                        paste("rep", params$rep, sep = "_"), sep = "")
data_location <- paste(simulation_dir,
                       paste("perc", params$perc, sep = "_"),
                       sep = "/")

mkdir_info <- paste(data_location, "output/vts/alignments", sep = "/")

where2save <- paste(mkdir_info, ".", sep = "/")

write.table(mkdir_info, file = "mkdir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)

write.table(where2save, file = "dir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)


# You have to download IQ-TREE to run this script
# Change to the correct path of IQ-TREE on your computer
#Software <- "/Applications/iqtree-2.1.2-MacOSX/bin/iqtree2"
Software <- "iqtree"
maxCPU <- as.numeric(get_params[4])


#list files
pattern <- paste("_", seq_length, ".fasta", sep = "")
list_files <- list.files(paste(data_location, "output/vts/alignments", sep = "/"),
                         pattern = pattern, full.names = TRUE)

list_sampleTimes <- list.files(paste(data_location, "output/vts/W", sep = "/"),
                               pattern = ".RData", full.names = TRUE)

list_sampleTimes <- list.files("output/vts/W",
                               pattern = ".RData", full.names = TRUE)

#make dir to save iqtree
iqtree_dirname <- paste("output/vts/alignments/iqtree_results_", seq_length, sep = "")
if (!dir.exists(iqtree_dirname)) {
  dir.create(iqtree_dirname, recursive = TRUE)
}

for(ali in list_files){


  ali_name <- paste("-s", ali, sep = " ")
  # -czb Collapse near zero branches, so that the final tree may be
  # multifurcating.
  iqtree_param <- c(ali_name, "-m HKY", "-T AUTO", "-ntmax", maxCPU, "-czb")
  #run iqtree
  system2(command = Software, args = iqtree_param)

  # get sampleTimes, cd4s, ehis to estimate treedater and infector probabilities
  load(list_sampleTimes)

  #estimate treedated tree
  #read ML tree
  mltree <- read.tree(paste(ali, ".treefile", sep = ""))

  #check if tree has polytomies
  #if tree has polytomies resolve polytomies randomly
  if(is.binary(mltree) == FALSE){
    mltree <- multi2di(mltree)
    if(is.rooted(mltree) == TRUE){
      mltree <- unroot(mltree)
    }
  }


  #run treedater
  dated_tree <- dater(tre = mltree, sts = sampleTimes[mltree$tip.label],
                      s = seqlength, clock = "uncorrelated")


  #drop tips of the tree that are from "global"
  # to calculate infector probability
  migrant_ID <- unlist(lapply(dated_tree$tip.label, function(x)
    ifelse(str_split(x, "_")[[1]][2] == 1 | str_split(x, "_")[[1]][2] == 21, FALSE, TRUE)))
  migrant_ID <- setNames(migrant_ID, dated_tree$tip.label)
  toDrop <- migrant_ID[migrant_ID == TRUE]
  region_only_dated_tree <- drop.tip(dated_tree, names(toDrop))

  # calculate infector probability (W) on the estimated tree
  # using same cd4s, ehis and MH as in the true trees
  # sampleTimes: must use years
  # cd4s: named numeric vector, cd4 at time of sampling
  # ehi: named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  # numberPeopleLivingWithHIV: scalar
  # numberNewInfectionsPerYear: scalar

  W_estimated <- phylo.source.attribution.hiv.msm( region_only_dated_tree, sampleTimes[region_only_dated_tree$tip.label],
                                                   cd4s = all_cd4s[region_only_dated_tree$tip.label],
                                                   ehi = ehis[region_only_dated_tree$tip.label],
                                                   numberPeopleLivingWithHIV  = totalPLWHIV,
                                                   numberNewInfectionsPerYear = newinf_per_year,
                                                   maxHeight = years,
                                                   res = 1e3,
                                                   treeErrorTol = Inf)

  #Create directory named W_estimated (to save everything related to infector probability)
  # if it does not exist
  if (!dir.exists("output/vts/W_estimated")) {
    dir.create("output/vts/W_estimated")
  }

  W_filename <- paste("output/vts/W_estimated/", seq_length,
                      "_migrant_years_1_simple_estimated", ".RData", sep="")

  save(years, last_sample_date, tm, dated_tree, region_only_dated_tree,
       sampleTimes, all_cd4s, ehis, newinf_per_year, totalPLWHIV, W_estimated,
       file = W_filename)
}

#move iqtree results to directory for iqtree
#files <- paste(getwd(), "output/vts/alignments", "*.fasta.*", sep = "/")
files <- paste(mkdir_info, "*.fasta.*", sep = "/")
#iqtree_dir <- paste(getwd(), iqtree_dirname, sep = "/")
iqtree_dir <- paste(data_location, iqtree_dirname, sep = "/")

#mkdir in ephemeral
if (!dir.exists(iqtree_dir)) {
  dir.create(iqtree_dir, recursive = TRUE)
}
#command_args <- paste("mv", files, iqtree_dir, sep = " ")
command_args <- paste("mv", files, paste(iqtree_dir, ".", sep = "/"), sep = " ")
system(command_args)


summaryW(sim = paste(params$params, params$rep, params$perc, sep = "."),
         tm = tm, W_estimated,
         tree = region_only_dated_tree, code = seq_length,
         prefix = paste("output/vts/W_estimated", seq_length, sep = "/"),
         labels = TRUE)

#end of script
end_time <- Sys.time()
print("IQTREE simulation took:")
end_time - start_time

iqtree_time <- data.frame(start = start_time, end = end_time)
saveRDS(iqtree_time, paste(seq_length, "iqtree_time.RDS", sep = "_"))
