# create pbs script to run shiver in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

#number of subjobs in array job
n = 10000

arguments <- commandArgs(trailingOnly = TRUE)

#line_number <- as.numeric(arguments[1])
#line_number <- 1
#params <- readRDS(system.file("data/simulations_deepseq_cluster.RDS",
#                              package = "HIVepisimAnalysis"))[line_number,]

#if in Imperial College cluster
cluster <- "/rds/general/user/fferre15/ephemeral/deepseq"
#if in XSEDE cluster
#cluster <- "/ocean/projects/bio210082p/nascimen/data_deepsequencing/"

#cluster <- as.character(arguments[1])

migration_variation <- as.character(arguments[1])
#migration_variation <- "50migrants"
migration_variation <- paste("best_trajectories", migration_variation, sep ="_")



simulation_dir <- paste(cluster, migration_variation, sep = "/")

simulation_dirs <- dir(dir(simulation_dir,full.names=TRUE), full.names=TRUE)

#location of sequence alignments
fasta_loc <- "seqgen/alignments/by_ID/results_sampling/amplicons"

outfiles <- list.files(paste(simulation_dirs, fasta_loc, sep = "/"), full.names = TRUE, pattern = "fasta")

input_split_dirs <- split_dirs(outfiles, size = n)

output_files <- str_split(outfiles, pattern = "/")
output_files <- unlist(lapply(output_files, function(x) paste(x[1:10], collapse = "/")))
output_split_dirs <- split_dirs(output_files, size = n)

#mknew dir to copy Illumina read files
mkdir_outfiles <- paste(output_files, "Illumina_reads", sep = "/")
mkdir_split_dirs <- split_dirs(mkdir_outfiles, size = n)



#location of alignments by ID
#ali_by_ID_loc <- "output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling"



#where2save <- paste(mkdir_info, ".", sep = "/")

#write.table(mkdir_info, file = "mkdir_name.txt",
#            quote = FALSE, row.names = FALSE,
#            col.names = FALSE)

#write.table(where2save, file = "dir_name.txt",
#            quote = FALSE, row.names = FALSE,
#            col.names = FALSE)


for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list_artIllumina", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  mkdir_filename <- paste(paste("mkdir_file_list_artIllumina", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = mkdir_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)


  pbsfilename <- paste(paste("ART", i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_dirs[[i]]), sep = "")


  pbstext <- paste("#PBS -l walltime=8:00:00",
                   "#PBS -l select=1:ncpus=1:mem=10gb",
                   paste("#PBS -o ", "Illumina", "_ART.stdout", sep = ""),
                   paste("#PBS -e ", "Illumina", "_ART.stderr", sep = ""),
                   array_number,
                   sep = "\n")

  pbstext <- paste(pbstext,
                   "\n",
                   "## load in the R environment",
                   "module load anaconda3/personal",
                   "source activate deep_seq_analysis",

                   "\n",
                   "## copy required files to the temporary directory on the compute node",
                   "export JOB_NUM=$(echo ${PBS_JOBID} | cut -f 1 -d '.' | cut -f 1 -d '[')",
                   "export WORKDIR=\"${EPHEMERAL}/${JOB_NUM}.${PBS_ARRAY_INDEX}\"",
                   "mkdir -p $WORKDIR",
                   "cd $WORKDIR",

                   "\n",
                   "#copy R scripts",
                   "cp $HOME/Ethics-HIV/large_pop/Rscripts/deepsequencing/3.Simulate_Illumina_reads_v4.R .",

                   "\n",

                   "cp $HOME/Ethics-HIV/Programs/artbinmountrainier2016.06.05linux64.tgz .",

                   "\n",
                   "#untar ART_Illumina program",
                   "tar -xzvf artbinmountrainier2016.06.05linux64.tgz ",
                   sep = "\n")

  #copy files to run shiver
  #create dirs to copy files
  new_art_dir <- "Illumina_reads"
  #new_fasta_dir <- "output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling/amplicons"


  pbstext <- paste(pbstext,


                   "#make new dir",
                   paste("mkdir", "-p", fasta_loc, sep = " "),

                   "\n",
                   "#cp fasta sequences to run ART",
                   paste(paste("cp -a ", "$(head -n $PBS_ARRAY_INDEX ", input_filename, "| tail -1)", sep = ""),
                         paste(fasta_loc, "/.", sep = "")),
                   "\n",
                   "## Run R scripts",
                   "Rscript 3.Simulate_Illumina_reads_v4.R",


                   "\n",
                   "#Remove directories for ART",
                   "# so they do not get transferred back to home directory",
                   "rm -rf art_bin_MountRainier",
                   "rm artbinmountrainier2016.06.05linux64.tgz",

                   "\n",
                   "## make dir in the submission directory - anything not copied back will be lost",
                   paste("mkdir -p", "$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = " "),

                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", paste("$WORKDIR", new_art_dir, "* ", sep = "/"), "$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""),

                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}

