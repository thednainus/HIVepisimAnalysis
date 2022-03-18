# create job script to run ART Illumina in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

arguments <- commandArgs(trailingOnly = TRUE)

#line_number <- as.numeric(arguments[1])
#line_number <- 1
#params <- readRDS(system.file("data/simulations_deepseq_cluster.RDS",
#                              package = "HIVepisimAnalysis"))[line_number,]

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


#location of alignments by ID
ali_by_ID_loc <- "output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling"


#copy files to run shiver
#create dirs to copy files
new_art_dir <- "Illumina_reads/results_sampling"
new_fasta_dir <- "output_deepseq/vts/merged_trees/alignments/by_ID"


#input files to run ART Illumina
input_files <- dir(paste(simulation_dirs, ali_by_ID_loc, sep = "/"), full.names=TRUE)

output_files <- str_split(input_files, pattern = "/")
output_files <- unlist(lapply(output_files, function(x) paste(x[1:10], collapse = "/")))


#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs
input_split_files <- split_dirs(input_files, size = 10000)

mkdir_split_files <- split_dirs(paste(output_files, new_art_dir, sep = "/"), size = 10000)



for (i in 1:length(input_split_files)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list_artIllumina", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_files[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  mkdir_filename <- paste(paste("mkdir_file_list_artIllumina", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = mkdir_split_files[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)


  jobfilename <- paste(paste("ART", i, sep = "_"), ".job", sep = "")

  array_number <- paste("#SBATCH", " --array=", "1-",
                          length(input_split_files[[i]]), sep = "")

  array_number <- paste("#SBATCH", " --array=", "1-",
                        length(input_split_files[[i]]), sep = "")

  jobtext <- paste("#!/bin/bash",
                   "#SBATCH -t 01:00:00",
                   "#SBATCH -p RM-shared",
                   "#SBATCH --ntasks-per-node 1",
                   "#SBATCH --mem=2000",
                   array_number,
                   sep = "\n")



  jobtext <- paste(jobtext,
                   "\n",
                   "#echo commands to stdout",
                   "set -x",
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
                   "cp $HOME/large_pop/Rscripts/3.Simulate_Illumina_reads.R .",
                   paste("cp $HOME/large_pop/deep_sequencing/", input_filename, sep = ""),
                   paste("cp $HOME/large_pop/deep_sequencing/", mkdir_filename, sep = ""),


                   "\n",
                   sep = "\n")


  jobtext <- paste(jobtext,


                   "#make new dir",
                   paste("mkdir", "-p", new_fasta_dir, sep = " "),

                   "\n",
                   "#cp fasta sequences to run ART",
                   paste(paste("cp -a ", "$(head -n $SLURM_ARRAY_TASK_ID ", input_filename, "| tail -1)", sep = ""),
                         paste(new_fasta_dir, "/.", sep = "")),
                   "\n",
                   "## Run R scripts",
                   "Rscript 3.Simulate_Illumina_reads.R",


                   "\n",
                   "##make dir in $PROJECT to copy output files from ART Illumina",
                   paste("mkdir -p", "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)", sep = " "),
                   "\n",
                   "## copy files back to $PROJECT directory",
                   paste(paste("cp -a ", paste("$WORKDIR", ali_by_ID_loc, "Illumina_reads", "* ", sep = "/"),
                               "$(head -n $SLURM_ARRAY_TASK_ID ", mkdir_filename, "| tail -1)", sep = ""),
                         "/.", sep = ""),

                   sep = "\n")


  write.table(jobtext, file = jobfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}

