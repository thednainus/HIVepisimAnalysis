# prepare files to run phyloscanner

# shiver output files were run in individual directories per ID (shiver should
# be run in individual directory per ID), but in the cluster.

# for the cluster analysis, I moved the results from shiver to a directory
#named "../phyloscanner"


#start of script
start_time <- Sys.time()

library(DescTools)
library(stringr)

#prefix for path names
prefix <- paste(getwd(), "/", sep = "")

#phyloscanner_dir <- paste(prefix, "output_deepseq/vts/merged_trees/Illumina_reads/results_merged_trees_active_region_diag/phyloscanner",
#                         sep = "")
#phyloscanner_dir <- paste(prefix, "output_deepseq/vts/merged_trees/Illumina_reads/results_merged_trees_active_region_diag/phyloscanner",
#                         sep = "")
phyloscanner_dir <- paste(prefix, "phyloscanner_results", sep = "")

# Create input csv file for phyloscanner ----
#list remap_bam files
bam_files = list.files(phyloscanner_dir, pattern="*._remap.bam$")
ref_files = list.files(phyloscanner_dir, pattern="*.remap_ref.fasta")


if(length(bam_files) != length(ref_files)){
  stop("length of bam_files has to be equal to length of ref_files")
}

# create dataframe of bam and ref files

input_files <- data.frame(bam = bam_files, ref_fasta = ref_files)

# check that ID names in bam matches ID names in ref_fasta

input_files["bam_ids"] <- unlist(lapply(input_files$bam, function(x) str_split(x, pattern = "_")[[1]][2]))
input_files["ref_ids"] <- unlist(lapply(input_files$ref_fasta, function(x) str_split(x, pattern = "_")[[1]][2]))

# check that all IDs in order are equal to each other
if(any(input_files$bam_ids == input_files$ref_ids) == FALSE){
  stop("All or some IDs between bam and ref.fasta files do not match")
}

input_files["IDs"] <- paste("ID", input_files$bam_ids, sep = "_")

#paste full file name path to bam and ref.fasta files
input_files["bam"] <- paste(phyloscanner_dir, input_files$bam, sep = "/")
input_files["ref_fasta"] <- paste(phyloscanner_dir, input_files$ref_fasta, sep = "/")

#final dataframe
input_files <- input_files[c(1,2,5)]

#save csv file that will be used as input file for phyloscanner
filename <- paste(phyloscanner_dir, "BamsRefsAndIDs.csv", sep = "/")
write.table(x =  input_files, file = filename, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = FALSE)


#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time


prepare_phyloscanner <- data.frame(start = start_time, end = end_time)
saveRDS(prepare_phyloscanner, "prepare_phyloscanner.RDS")
