# create pbs script to run shiver in individual IDs.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

#  commandArgs(trailingOnly = TRUE) allow to call Rscript with an argument
# this argument will be in the form of "sim1" or "sim2", etc.
simulation <- commandArgs(trailingOnly = TRUE)

#outfiles <- list.files(paste("/rds/general/user/fferre15/ephemeral/res_", simulation, "/stage2", sep = ""), full.names = TRUE)
#outfiles <- list.files("/rds/general/user/fferre15/home/Ethics-HIV/small_pop/sim1/res_sim1/stage2", full.names = TRUE)
#outfiles <- list.files("output_deepseq1/vts/merged_trees/Illumina_reads/results_sampling/")
outfiles <- list.files("/Users/user/Desktop/teste/storage_test", full.names = TRUE)



#save input list file
input_files <- list()

for(i in 1:length(outfiles)){

  dir_list <- list.dirs(outfiles[i])
  #dir_list <- dir_list[grepl("ID_", dir_list)]
  dir_list <- dir_list[grepl("/[0:9]", dir_list)]
  input_files[[i]] <- paste(dir_list, "/.", sep = "")

}

input_files <- unlist(input_files)

#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(input_files, size = 10000)


output_files <- str_split(input_files, pattern = "/")
output_files <- unlist(lapply(output_files, function(x) paste(x[1:7], collapse = "/")))
#output_files <- unlist(lapply(output_files, function(x) paste(x[1:9], collapse = "/")))
output_split_dirs <- split_dirs(output_files, size = 10000)

mkdir_outfiles <- lapply(output_split_dirs, function(x) str_split(x, pattern = "/"))

mkdir_outfiles <- lapply(mkdir_outfiles, lapply, function(x) str_split(x, pattern = "/")[[7]])
#mkdir_outfiles <- lapply(mkdir_outfiles, lapply, function(x) str_split(x, pattern = "/")[[9]])
mkdir_outfiles <- lapply(mkdir_outfiles, function(x) unlist(x))


for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list", simulation, i, sep ="_"),
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
  mkdir_filename <- paste(paste("mkdir_file_list", simulation, i, sep ="_"),
                           "txt", sep = ".")
  write.table(x = mkdir_outfiles[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)

  mknew_dir <- paste("$(head -n $PBS_ARRAY_INDEX", mkdir_filename, "| tail -1)", sep = " ")


  pbsfilename <- paste(paste("shiver", simulation, i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_dirs[[i]]), sep = "")


  pbstext <- paste("#PBS -l walltime=72:00:00",
                   "#PBS -l select=1:ncpus=1:mem=10gb",
                   paste("#PBS -o ", simulation, "_stage3.stdout", sep = ""),
                   paste("#PBS -e ", simulation, "_stage3.stderr", sep = ""),
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
                   "cp $HOME/Ethics-HIV/small_pop/Rscripts/4.Map_reads_shiver_v2.R .",

                   "\n",
                   "#copy Softwares and other input files",
                   "#that we will need to run the abovementioned R scripts",
                   "cp $HOME/Ethics-HIV/small_pop/input_files/reference_subtypeB_completeGenome.fasta .",

                   "cp $HOME/Ethics-HIV/Programs/shiver.tar.gz .",

                   paste("cp $HOME/Ethics-HIV/small_pop/", simulation, "/", input_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/small_pop/", simulation, "/", mkdir_filename, " .", sep = ""),

                   "\n",
                   "#untar shiver program",
                   "tar -xzvf shiver.tar.gz ",
                   sep = "\n")

  #copy files to run shiver
  #create dirs to copy files
  new_shiver_dir <- "shiver2"
  new_reads_dir <- paste(new_shiver_dir, "Illumina_reads", sep = "/")
  output_dirname <- "phyloscanner"

  new_ephemeraldir <- paste("$EPHEMERAL",paste("shiver", simulation, sep = "_"),
                            mknew_dir, output_dirname,  sep = "/")

  pbstext <- paste(pbstext,


                   "#make new dir",
                   paste("mkdir", "-p", new_reads_dir, sep = " "),

                   "\n",
                   "#cp Illumina reads to run shiver",
                   paste(paste("cp -a ", "$(head -n $PBS_ARRAY_INDEX ", input_filename, "| tail -1)", sep = ""),
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
                   "#Remove directories for shiver",
                   "# so they do not get transferred back to home directory",
                   "rm -rf shiver",
                   "rm shiver.tar.gz",

                   "\n",
                   "## make dir in the submission directory - anything not copied back will be lost",
                   paste("mkdir -p", new_ephemeraldir, sep = " "),

                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", "$WORKDIR/shiver2/*.bam ", new_ephemeraldir, "/.", sep = ""),
                   paste("cp -a ", "$WORKDIR/shiver2/*.bai ", new_ephemeraldir, "/.", sep = ""),
                   paste("cp -a ", "$WORKDIR/shiver2/*ref.fasta ", new_ephemeraldir, "/.", sep = ""),

                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}
























