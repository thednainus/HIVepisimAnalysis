# script to run phyloscanner: https://github.com/BDI-pathogens/phyloscanner
# this script run the part of phyloscanner to generate read alignments:
# we use the command in phyloscanner --no-trees to supress the reconstruction
# of phylogenetic trees

# phyloscanner dependencies
# https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh

#start of script
start_time <- Sys.time()

library(DescTools)
library(reticulate)

use_python("/Users/user/opt/miniconda2/bin/python2")
#use_python("/rds/general/user/fferre15/home/anaconda3/envs/deep_seq_analysis/bin/python2")
py_config()

#location of phyloscanner make tree python code
phyloscanner <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner"
#phyloscanner <- paste(prefix, "phyloscanner/", sep = "")
make_trees <- paste(phyloscanner, "phyloscanner_make_trees.py", sep = "/")

# window parameters
# these window parameters are from Table S1 from Ratman et al. 2019 (Nature Communications)
#window_width <- "250,"
#window_overlap <- "125,"
#start_pos <- "800,"
#end_pos <- "9400"
#get argument of window to generate alignmets
#in the form of 800,1049 (will generate alignment from 800 to 1049)
#read csv with window coordinate
row_number <- commandArgs(trailingOnly = TRUE)
row_number <- as.numeric(row_number)
#row_number <- 1
windows <- read.table(system.file("data/phyloscanner_windows_30Sep2021.csv",
                                package = "HIVepisimAnalysis"))[row_number,]
windows <- paste("--windows", windows, sep = " ")



# make read alignments
#alignment_other_refs <- paste("--alignment-of-other-refs",
#                              paste(phyloscanner, "InfoAndInputs",
#                                    "2refs_HXB2_C.BW.fasta", sep = "/"), sep = " ")

alignment_other_refs <- paste("--alignment-of-other-refs",
                              "HIV1_REF_2020_phyloscanner.fasta", sep = " ")

paiwise_align_to <- paste("--pairwise-align-to",
                          "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                          sep = " ")

merge_pairs <- "--merge-paired-reads"

quality_trim_ends <- paste("--quality-trim-ends", "23", sep = " ")

min_internal_quality <- paste("--min-internal-quality", "23", sep = " ")

excision_ref <- paste("--excision-ref", "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                      sep = " " )

excision_coords <- paste("--excision-coords", "$(cat",
                         paste(phyloscanner, "InfoAndInputs",
                               "DrugResistancePositionsInHXB2.txt", sep = "/"),
                         ")", sep = " ")


merging_threshold <- paste("--merging-threshold-a", "1", sep = " ")

min_read_count <- paste("--min-read-count", "2", sep = " ")

#option to generate read alignments only
#trees will not be generated
no_trees <- "--no-trees"


#Phyloscanner location

# list directories
# phyloscanner directory (with the bam files)
phyloscanner_dir_fullPath <- paste(getwd(), "phyloscanner", sep = "/")

BamsRefsAndIds <- paste(phyloscanner_dir_fullPath, "BamsRefsAndIDs.csv", sep = "/")


#parameters <- paste(BamsRefsAndIds, windows, alignment_other_refs,
#              paiwise_align_to, merge_pairs, quality_trim_ends,
#              min_internal_quality, excision_ref, excision_coords,
#              merging_threshold, min_read_count, raxml, sep = " ")

parameters <- paste(BamsRefsAndIds, windows, alignment_other_refs,
                    paiwise_align_to, merge_pairs, quality_trim_ends,
                    min_internal_quality,
                    excision_ref, excision_coords,
                    merging_threshold, min_read_count, no_trees, sep = " ")


command_make_trees <- paste("cd", phyloscanner_dir_fullPath, "&&", make_trees, sep = " ")
makeTrees_and_args <- paste(command_make_trees, parameters, sep = " ")


system(makeTrees_and_args)



#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time


read_alignments_time <- data.frame(start = start_time, end = end_time)
saveRDS(read_alignments_time, paste(as.character(row_number),
                                    "phyloscanner_read_alignments_time.RDS",
                                    sep = "_"))

