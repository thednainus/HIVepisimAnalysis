# create pbs script to run iqtree analysis on fasta files for reads
# to analyse with phyloscanner.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

dir_name <- "/Users/user/Desktop/teste/phyloscanner_results"

phyloscanner_results <- list.files(dir_name, full.names = TRUE)

#save input list file
input_files <- list()

for(i in 1:length(phyloscanner_results)){


  output_files <- paste(phyloscanner_results[i], "phyloscanner_results/AlignedReads", sep = "/")

  outfiles_excised <- list.files(output_files, pattern = "PositionsExcised")
  outfiles <- list.files(output_files, pattern = "AlignedReadsInWindow_\\d+_to_\\d+")


  # get filenames to estimate trees
  # if a filename that contains the word "PositionsExcised" exists, then I keep
  # this file and ignore the other file from the same

  toDelete <- str_split(outfiles_excised, pattern = "_")
  toDelete <- unlist(lapply(toDelete, function(x) paste(x[c(1,3,4,5)], collapse = "_")))

  toKeep <- outfiles[! outfiles %in% toDelete]

  #merge filenames to estimate phylogenetic trees
  all_files <- c(toKeep, outfiles_excised)
  input_files[[i]] <- paste(phyloscanner_results[i], all_files, sep = "/")

}


input_files <- unlist(input_files)

#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(input_files, size = 10000)


output_files <- str_split(input_files, pattern = "/")
output_files <- unlist(lapply(output_files, function(x) paste(x[1:6], collapse = "/")))
#output_files <- unlist(lapply(output_files, function(x) paste(x[1:8], collapse = "/")))
output_split_dirs <- split_dirs(output_files, size = 10000)

mkdir_outfiles <- lapply(output_split_dirs, function(x) str_split(x, pattern = "/"))
mkdir_outfiles <- lapply(mkdir_outfiles, lapply, function(x) str_split(x, pattern = "/")[[5]])
#mkdir_outfiles <- lapply(mkdir_outfiles, lapply, function(x) str_split(x, pattern = "/")[[8]])
mkdir_outfiles <- lapply(mkdir_outfiles, function(x) unlist(x))




simulation <- "sim1"

#save input list file

#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(input_files, size = 10000)


for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_files_iqtree", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  #create files with mkdir file list
  mkdir_filename <- paste(paste("mkdir_files_iqtree", simulation, i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = mkdir_outfiles[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)



  mknew_dir <- paste("$(head -n $PBS_ARRAY_INDEX", mkdir_filename, "| tail -1)", sep = " ")


  pbsfilename <- paste(paste("iqtree_reads", simulation, i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_dirs[[i]]), sep = "")




  pbstext <- paste("#PBS -l walltime=72:00:00",
                   "#PBS -l select=1:ncpus=8:mem=10gb",
                   "#PBS -o sim1_iqtree.stdout",
                   "#PBS -e sim1_iqtree.stderr",
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
                   "cp $HOME/Ethics-HIV/small_pop/Rscripts/7.Estimate_trees_reads_IQTREE.R .",

                   paste("cp $HOME/Ethics-HIV/small_pop/sim1/", input_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/small_pop/sim1/", mkdir_filename, " .", sep = ""),

                  sep = "\n")

  new_ephemeraldir <- paste("$EPHEMERAL/shiver_sim1", mknew_dir, "phyloscanner_results/iqtree/.", sep = "/")

  pbstext <- paste(pbstext,

                   "\n",
                   "#cp aligned read to run IQTREE",
                   paste(paste("cp -a ", "$(head -n $PBS_ARRAY_INDEX ", input_filename, "| tail -1)", sep = ""), ".", sep = " "),

                   "\n",
                   "## Run R scripts",
                   "Rscript 7.Estimate_trees_reads_IQTREE.R",

                   "\n",
                   "# remove fasta file",
                   paste("rm ", "*.fasta", sep = ""),

                   "\n",
                   "## make dir in the submission directory - anything not copied back will be lost",
                   paste("mkdir -p", new_ephemeraldir, sep = " "),


                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", "$WORKDIR/.", new_ephemeraldir, sep = " "),

                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}
