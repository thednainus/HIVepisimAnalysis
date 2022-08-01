# create pbs script to run shiver in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

arguments <- commandArgs(trailingOnly = TRUE)

line_number <- as.numeric(arguments[1])
#line_number <- 1
params <- readRDS(system.file("data/simulations_deepseq_cluster.RDS",
                              package = "HIVepisimAnalysis"))[line_number,]

#if in Imperial College cluster
#cluster <- "/rds/general/user/fferre15/ephemeral/"
#if in XSEDE cluster
#cluster <- "/ocean/projects/bio210082p/nascimen/data_deepsequencing/"

cluster <- as.character(arguments[2])

migration_variation <- as.character(arguments[3])
#migration_variation <- "50migrants"
migration_variation <- paste("best_trajectories", migration_variation, sep ="_")



simulation_dir <- paste(paste(cluster, migration_variation, "/", sep = ""),
                        paste("params", params$params, sep = "_"),
                        "/",
                        paste("rep", params$rep, sep = "_"), sep = "")

mkdir_info <- paste(simulation_dir,
                    paste("perc", params$perc, sep = "_"),
                    sep = "/")

#location of alignments by ID
ali_by_ID_loc <- "output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling"



where2save <- paste(mkdir_info, ".", sep = "/")

write.table(mkdir_info, file = "mkdir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)

write.table(where2save, file = "dir_name.txt",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)


#outfiles <- list.files(paste("/rds/general/user/fferre15/ephemeral/res_", simulation, "/stage2", sep = ""), full.names = TRUE)
#outfiles <- list.files("/rds/general/user/fferre15/home/Ethics-HIV/small_pop/sim1/res_sim1/stage2", full.names = TRUE)
#outfiles <- list.files("output_deepseq1/vts/merged_trees/Illumina_reads/results_sampling/")

#input files to run ART Illumina
input_files <- list.files(paste(mkdir_info, ali_by_ID_loc, sep = "/"), full.names = TRUE)
#input_files <- list.files("/Users/user/Desktop/teste/ART_teste/output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling", full.names = TRUE)


#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs
input_split_files <- split_dirs(input_files, size = 1000)



for (i in 1:length(input_split_files)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list_artIllumina", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_files[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)


  pbsfilename <- paste(paste("ART", simulation, i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_files[[i]]), sep = "")


  pbstext <- paste("#PBS -l walltime=72:00:00",
                   "#PBS -l select=1:ncpus=1:mem=10gb",
                   paste("#PBS -o ", simulation, "_ART.stdout", sep = ""),
                   paste("#PBS -e ", simulation, "_ART.stderr", sep = ""),
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
                   "cp $HOME/Ethics-HIV/small_pop/Rscripts/3.Simulate_Illumina_reads.R .",

                   "\n",

                   "cp $HOME/Ethics-HIV/Programs/artbinmountrainier2016.06.05linux64.tgz .",

                   paste("cp $HOME/Ethics-HIV/small_pop/deep_sequencing/", simulation, "/", input_filename, " .", sep = ""),
                   #paste("cp $HOME/Ethics-HIV/small_pop/deep_sequencing/", simulation, "/", mkdir_filename, " .", sep = ""),

                   "\n",
                   "#untar shiver program",
                   "tar -xzvf artbinmountrainier2016.06.05linux64.tgz ",
                   sep = "\n")

  #copy files to run shiver
  #create dirs to copy files
  new_art_dir <- "Illumina_reads/results_sampling"
  new_fasta_dir <- "output_deepseq/vts/merged_trees/alignments/by_ID"


  new_ephemeraldir <- paste("$EPHEMERAL", paste("SanDiego_deepseq/stage2/output_deepseq/vts/merged_trees", sep = "_"),
                            new_art_dir,  sep = "/")

  pbstext <- paste(pbstext,


                   "#make new dir",
                   paste("mkdir", "-p", new_fasta_dir, sep = " "),

                   "\n",
                   "#cp fasta sequences to run ART",
                   paste(paste("cp -a ", "$(head -n $PBS_ARRAY_INDEX ", input_filename, "| tail -1)", sep = ""),
                         paste(new_fasta_dir, "/.", sep = "")),
                   "\n",
                   "## Run R scripts",
                   "Rscript 3.Simulate_Illumina_reads.R",


                   "\n",
                   "#Remove directories for ART",
                   "# so they do not get transferred back to home directory",
                   "rm -rf art_bin_MountRainier",
                   "rm artbinmountrainier2016.06.05linux64.tgz",

                   "\n",
                   "## make dir in the submission directory - anything not copied back will be lost",
                   paste("mkdir -p", new_ephemeraldir, sep = " "),

                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", paste("$WORKDIR", new_art_dir, "* ", sep = "/"), new_ephemeraldir, "/.", sep = ""),

                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}

