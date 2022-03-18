# create pbs script to run shiver in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

#number of subjobs in array job
n = 1000

#  commandArgs(trailingOnly = TRUE) allow to call Rscript with an argument
# this argument will be in the form of "sim1" or "sim2", etc.
arguments <- commandArgs(trailingOnly = TRUE)


#if in Imperial College cluster
#cluster <- "/rds/general/user/fferre15/ephemeral/"
#if in XSEDE cluster
#cluster <- "/ocean/projects/bio210082p/nascimen/data_deepsequencing/"

cluster <- as.character(arguments[1])

migration_variation <- as.character(arguments[2])
#migration_variation <- "50migrants"
migration_variation <- paste("best_trajectories", migration_variation, sep ="_")

simulation_dirs <- dir(dir(dir(paste(cluster, migration_variation, sep = ""),
                               full.names=TRUE), full.names=TRUE),
                       pattern="perc*", full.names=TRUE)


#location of ART Illumina results
ART_Illumina_loc <- "Illumina_reads/results_sampling"

#output files from Illumina ART simulations
#outfiles <- list.files(paste("/rds/general/user/fferre15/ephemeral/res_", simulation, "/stage2", sep = ""), full.names = TRUE)
#outfiles <- list.files("/rds/general/user/fferre15/home/Ethics-HIV/small_pop/sim1/res_sim1/stage2", full.names = TRUE)
#outfiles <- list.files("output_deepseq1/vts/merged_trees/Illumina_reads/results_sampling/")
outfiles <- list.files(paste(simulation_dirs, ART_Illumina_loc, sep = "/"), full.names = TRUE)



#save input list file
#so I have something like
#" "/ocean/projects/bio210082p/nascimen/data_deepsequencing/best_trajectories_50migrants/params_1067/rep_99/perc_0.3/Illumina_reads/results_sampling/ID_99981/."
input_files <- paste(outfiles, "/.", sep = "")

#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(input_files, size = n)


output_files <- str_split(input_files, pattern = "/")
#output_files <- unlist(lapply(output_files, function(x) paste(x[1:7], collapse = "/")))
#output_files <- unlist(lapply(output_files, function(x) paste(x[1:9], collapse = "/")))
output_files <- unlist(lapply(output_files, function(x) paste(x[1:10], collapse = "/")))
output_split_dirs <- split_dirs(output_files, size = n)

#mknew dir to copy shiver files
#mkdir_outfiles <- lapply(output_split_dirs, function(x) str_split(x, pattern = "/"))
mkdir_outfiles <- lapply(output_split_dirs, function(x) paste(x, "shiver_results", sep = "/"))



for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  #create files with output file list
  #output_filename <- paste(paste("output_file_list", simulation, i, sep ="_"),
  #                         "txt", sep = ".")
  #write.table(x = output_split_dirs[[i]], quote = FALSE,
  #            row.names = FALSE, col.names = FALSE,
  #            file = output_filename)

  #create files with mkdir file list
  mkdir_filename <- paste(paste("mkdir_file_list", i, sep ="_"),
                           "txt", sep = ".")
  write.table(x = mkdir_outfiles[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)

  mknew_dir <- paste("$(head -n $SLURM_ARRAY_TASK_ID", mkdir_filename, "| tail -1)", sep = " ")


  jobfilename <- paste(paste("shiver", i, sep = "_"), ".job", sep = "")


  array_number <- paste("#SBATCH", " --array=", "1-", length(input_split_dirs[[i]]), sep = "")


  jobtext <- paste("#!/bin/bash",
                   "#SBATCH -t 72:00:00",
                   "#SBATCH -p RM-shared",
                   "#SBATCH --ntasks-per-node 1",
                   "#SBATCH --mem=2000",
                   array_number,
                   sep = "\n")


  jobtext <- paste(jobtext,
                   "\n",
                   "## load in the R environment",
                   "module load anaconda3/2020.11",
                   "source activate HIVepisim",

                   "\n",
                   "## copy required files to the temporary directory on the compute node",
                   "WORKDIR=${PROJECT}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}",
                   "mkdir -p $WORKDIR",
                   "cd $WORKDIR",

                   "\n",
                   "#copy R scripts",
                   "cp $HOME/large_pop/Rscripts/4.Map_reads_shiver_v2.R .",

                   "\n",
                   "#copy other input files",
                   "#that we will need to run the abovementioned R scripts",
                   "cp $HOME/large_pop/deep_sequencing/reference_subtypeB_completeGenome.fasta .",

                   paste("cp $HOME/large_pop/deep_sequencing/", input_filename, " .", sep = "/"),
                   paste("cp $HOME/large_pop/deep_sequencing/", mkdir_filename, " .", sep = "/"),

                   sep = "\n")

  #copy files to run shiver
  #create dirs to copy files
  new_shiver_dir <- "shiver2"
  new_reads_dir <- paste(new_shiver_dir, "Illumina_reads", sep = "/")
  output_dirname <- "phyloscanner"

  jobtext <- paste(jobtext,

                   "\n",

                   "#make new dir",
                   paste("mkdir", "-p", new_reads_dir, sep = " "),

                   "\n",
                   "#cp Illumina reads to run shiver",
                   paste(paste("cp -a ", "$(head -n $SLURM_ARRAY_TASK_ID ", input_filename, "| tail -1)", sep = ""),
                         paste(new_reads_dir, "/.", sep = "")),
                   "\n",
                   "## Run R scripts",
                   "Rscript 4.Map_reads_shiver_v2.R",

                   "\n",
                   "# remove simulated Illumina reads to decrease storage spage",
                   paste("rm ", new_shiver_dir, "/*.fq", sep = ""),
                   paste("rm ", new_reads_dir, "/*.fq.gz", sep = ""),

                   "\n",
                   "#remove tmp files generated by shiver",
                   paste("rm ", new_shiver_dir, "/temp_*", sep = ""),


                   "\n",
                   "#mkdir files to copy results of shiver to run phyloscanner later",
                   paste("mkdir -p", "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)", sep = " "),

                   "\n",

                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", "$WORKDIR/shiver2/*.bam ", "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)"),
                   paste("cp -a ", "$WORKDIR/shiver2/*.bai ", "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)"),
                   paste("cp -a ", "$WORKDIR/shiver2/*ref.fasta ", "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)"),

                   sep = "\n")


  write.table(jobtext, file = jobfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}
