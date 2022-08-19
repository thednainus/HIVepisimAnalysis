# script to run phyloscanner: https://github.com/BDI-pathogens/phyloscanner
# this script run the part of phyloscanner to generate read alignments:
# we use the command in phyloscanner --no-trees to suppress the reconstruction
# of phylogenetic trees

# phyloscanner dependencies
# https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh

#start of script
start_time <- Sys.time()

library(DescTools)

params <- commandArgs(trailingOnly = TRUE)

#location of phyloscanner make tree python code
#phyloscanner <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner"
phyloscanner <- paste(getwd(), "phyloscanner", sep = "/")
make_trees <- paste(phyloscanner, "phyloscanner_make_trees.py", sep = "/")

prefix <- as.character(params[1])

line_number <- as.numeric(params[2])
reps <- list.files(list.files(prefix, full.names=T), full.names=T)[line_number]


alignment_other_refs <- paste("--alignment-of-other-refs",
                              "HIV1_REF_2020_phyloscanner.fasta", sep = " ")

paiwise_align_to <- paste("--pairwise-align-to",
                          "B.FR.83.HXB2_LAI_IIIB_BRU.K03455",
                          sep = " ")

merge_pairs <- "--merge-paired-reads"

quality_trim_ends <- paste("--quality-trim-ends", "25", sep = " ")

min_internal_quality <- paste("--min-internal-quality", "25", sep = " ")

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
phyloscanner_dir_fullPath <- paste(reps, "shiver_results", sep = "/")

BamsRefsAndIds <- paste(phyloscanner_dir_fullPath, "BamsRefsAndIDs.csv", sep = "/")


#parameters <- paste(BamsRefsAndIds, windows, alignment_other_refs,
#              paiwise_align_to, merge_pairs, quality_trim_ends,
#              min_internal_quality, excision_ref, excision_coords,
#              merging_threshold, min_read_count, raxml, sep = " ")

# window parameters
# these window parameters are from Zhang et al. 2021 (Clinical Infectious Disease)
#get argument of window to generate alignments
#in the form of 800,1049 (will generate alignment from 800 to 1049)
#read csv with window coordinate
window_width <- "340,"
window_overlap <- "170,"
start_pos <- "800,"
end_pos <- "9400"

windows <- paste("--auto-window-params ", window_width, window_overlap, start_pos, end_pos, sep = "")

parameters <- paste(BamsRefsAndIds, windows, alignment_other_refs,
                    paiwise_align_to, merge_pairs, quality_trim_ends,
                    min_internal_quality,
                    excision_ref, excision_coords,
                    merging_threshold, min_read_count, no_trees, sep = " ")


makeTrees_and_args <- paste(make_trees, parameters, sep = " ")

system(makeTrees_and_args)


#where to save files to mknew dir
where2save_files <- paste(reps, "phyloscanner_ali_35_340", sep = "/")
write.table(where2save_files, file = "where2save_files.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

