# create pbs script to run shiver in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

#  commandArg s(trailingOnly = TRUE) allow to call Rscript with an argument
# this argument will be in the form of "sim1" or "sim2", etc.
#simulation <- commandArgs(trailingOnly = TRUE)
simulation="5perc"

#outfiles <- list.files(paste("/rds/general/user/fferre15/ephemeral/res_", simulation, "/stage2", sep = ""), full.names = TRUE)
#outfiles <- list.files("/rds/general/user/fferre15/home/Ethics-HIV/small_pop/sim1/res_sim1/stage2", full.names = TRUE)
#outfiles <- list.files("output_deepseq1/vts/merged_trees/Illumina_reads/results_sampling/")

#input files to run ART Illumina
input_files <- list.files("/Users/user/Desktop/teste/ART_teste/output_deepseq/vts/merged_trees/alignments/by_ID/results_sampling", full.names = TRUE)

#mkdir_outfiles <- rep(1:length(input_files))
#mkdir_split_outfiles <- split_dirs(mkdir_outfiles, size = 1000)


#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs
input_split_files <- split_dirs(input_files, size = 1000)



for (i in 1:length(input_split_files)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list_art", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_files[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)


  #create files with mkdir file list
  #mkdir_filename <- paste(paste("mkdir_file_list_art", simulation, i, sep ="_"),
  #                         "txt", sep = ".")
  #write.table(x = mkdir_split_outfiles[[i]], quote = FALSE,
  #            row.names = FALSE, col.names = FALSE,
  #            file = mkdir_filename)

  #mknew_dir <- paste("$(head -n $PBS_ARRAY_INDEX", mkdir_filename, "| tail -1)", sep = " ")

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

