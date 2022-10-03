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

#params <- readRDS(system.file("data/simulations_imperial_cluster.RDS",
#                              package = "HIVepisimAnalysis"))

#params <- params[c(1:30, 101:130, 201:230, 301:330, 401:430,
#                   501:530, 601:630, 701:730, 801:830, 901:930,
#                   1001:1030, 1101:1130, 1201:1230, 1301:1330, 1401:1430,
#                   1501:1530, 1601:1630, 1701:1730, 1801:1830, 1901:1930),]

params <- readRDS(as.character(get_params[5]))
params <- params[line_number,]

seqlength <- as.numeric(get_params[2])
seq_length <- paste(seqlength, "bp", sep = "")

#get where the data is located
migration_variation <- as.character(get_params[3])
#migration_variation <- paste("best_trajectories", migration_variation, sep ="_")
simulation_dir <- paste("/rds/general/user/fferre15/ephemeral/",
                        migration_variation, "/",
                        paste("params", params$param, sep = "_"),
                        "/",
                        paste("rep", params$rep, sep = "_"), sep = "")
data_location <- paste(simulation_dir,
                       "sampler2",
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

#cp fasta file to $WORKDIR
system(paste("cp", list_files[1], getwd(), sep = " "))

list_sampleTimes <- list.files(paste(data_location, "output/vts/W", sep = "/"),
                               pattern = ".RData", full.names = TRUE)

#make dir to save iqtree
iqtree_dirname <- paste("output/vts/alignments/iqtree_results_", seq_length, sep = "")
if (!dir.exists(iqtree_dirname)) {
  dir.create(iqtree_dirname, recursive = TRUE)
}


ali <- str_split(list_files, pattern = "/")[[1]][15]


ali_name <- paste("-s", ali, sep = " ")
# -czb Collapse near zero branches, so that the final tree may be
# multifurcating.
iqtree_param <- c(ali_name, "-m HKY", "-T AUTO", "-ntmax", maxCPU, "-czb")
#run iqtree
system2(command = Software, args = iqtree_param)

#move iqtree results to directory for iqtree
system(paste("mv", "*.fasta.*", iqtree_dirname, sep = " "))


#end of script
end_time <- Sys.time()
print("IQTREE simulation took:")
end_time - start_time

iqtree_time <- data.frame(start = start_time, end = end_time)
saveRDS(iqtree_time, paste(seq_length, "iqtree_time.RDS", sep = "_"))
