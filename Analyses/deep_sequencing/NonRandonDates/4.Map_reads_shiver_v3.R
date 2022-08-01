# map simulated Illumina reads using shiver (https://github.com/ChrisHIV/shiver)
# this script uses subtype reference sequece for mapping

# INFORMATION:
# Because this is an analysis in the cluster, I will analyse each ID with shiver
# in separate. The simulated Illumina reads will be located in a directory
# named shiver2.
# There is a directory named shiver in which the program shiver is located.

# Note that the simulated Illumina reads have to be in a separated directory
# and should not be in the same directory as the working directory for shiver
# otherwise a error will be generated.


#start of simulation
start_time <- Sys.time()


library(DescTools)
library(stringr)

#prefix for path names
prefix <- paste(getwd(), "/", sep = "")

# list directories
reads_output <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                  full.names = TRUE)
# This will return the fastq files located in directory shiver2. These files
# will be copied using the pbs file.
#reads_output <- dir(path = "shiver2/Illumina_reads",
#                 full.names = TRUE)
reads_output <- reads_output[2]

#location of shiver
#shiver <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/shiver/"
shiver <- paste(prefix, "shiver/", sep = "")


# location for shiver align contigs
shiver_align_contigs <- paste(shiver, "shiver_align_contigs.sh", sep = "")

# location for shiver map reads
shiver_map_reads <- paste(shiver, "shiver_map_reads.sh", sep = "")

#location for MyInitDir
MyInitDir <- paste(shiver, "MyInitDir", sep = "")

#location of config.sh file
config_file <- paste(shiver, "config_MyChanges.sh", sep = "")


#subtypeB reference sequence for mapping
#loc_subtypeB <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/reference_subtypeB/genome/reference_subtypeB_completeGenome.fasta"
loc_subtypeB <- paste(prefix, "reference_subtypeB_completeGenome.fasta", sep = "")

# shiver_output is the directory in which we will save all outputs from the
# shiver analysis.
shiver_output <- "shiver2"
shiver_output <- paste(prefix, shiver_output, sep = "")

#Create directory named shiver2 if it does not exist
if (!dir.exists(shiver_output)) {
  dir.create(shiver_output)
}



#run shiver_align_contigs
# We need the ID for the sampled that shiver will use as prefix
# We will get this ID from the fastq filename that is in the form of ID_427_1.fq.gz

#list fq.gz files
SID_shiver <- list.files(reads_output, full.names = TRUE)[1]
#SID_shiver <- reads_output[1]
#SID_shiver <- str_split(SID_shiver, "/")[[1]][3]
SID_shiver <- str_split(SID_shiver, "/")[[1]][5]
#SID_shiver <- str_split(SID_shiver, "_")[[1]][2]
#SID_shiver <- paste("ID", SID_shiver, sep = "_")

parameters1 <- paste(MyInitDir, config_file, loc_subtypeB, SID_shiver, sep = " ")
command_shiver <- paste("cd", shiver_output, "&&", "bash", shiver_align_contigs, sep = " ")
shiver_and_args <- paste(command_shiver, parameters1, sep = " ")

system(shiver_and_args)

#run shiver_map_reads.sh
#check if the better alignment file exists
SID_wRefs <- paste(SID_shiver, "cut_wRefs.fasta", sep = "_")
SID_wRefs_fullPath <- paste(shiver_output, SID_wRefs, sep = "/")

if(file.exists(SID_wRefs_fullPath) == FALSE){
  SID_wRefs <- paste(SID_shiver, "raw_wRefs.fasta", sep = "_")
  SID_wRefs_fullPath <- paste(shiver_output, SID_wRefs, sep = "/")
}

SID_blast <- paste(SID_shiver, "blast", sep = ".")
SID_blast_fullPath <- paste(shiver_output, SID_blast, sep = "/")

reads_1 <- paste(SID_shiver, "1.fq.gz", sep = "_")
reads_1_fullPath <- paste(getwd(),reads_output, reads_1, sep = "/")

reads_2 <- paste(SID_shiver, "2.fq.gz", sep = "_")
#reads_2_fullPath <- paste(shiver_output, "Illumina_reads", reads_2, sep = "/")
reads_2_fullPath <- paste(getwd(), reads_output, reads_2, sep = "/")

parameters2 <- paste(MyInitDir, config_file, loc_subtypeB, SID_shiver,
                    SID_blast_fullPath, SID_wRefs_fullPath,
                    reads_1_fullPath, reads_2_fullPath,  sep = " ")

command_shiver_map <- paste("cd", shiver_output, "&&", "bash", shiver_map_reads, sep = " ")
shiver_and_args2 <- paste(command_shiver_map, parameters2, sep = " ")

system(shiver_and_args2)


#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time

shiver_time <- data.frame(start = start_time, end = end_time)
saveRDS(shiver_time, "shiver_time.RDS")
