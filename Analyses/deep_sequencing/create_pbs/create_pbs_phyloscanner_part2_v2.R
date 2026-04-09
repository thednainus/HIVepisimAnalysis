#This script will generate a pbs script for individuals trees
#for phyloscanner analyses
#When trying to analyse all trees at the same time, job did not complete in
#72 hours

library(HIVepisimAnalysis)
library(stringr)


iqtree_results <- list.files(list.files(list.files("/rds/general/user/fferre15/ephemeral/deepseq/best_trajectories_250migrants", full.names=T), full.names=T), full.names=T, pattern="iqtree")


input_files <- list()

for(i in 1:length(iqtree_results)){
  iqtree_results1 <- iqtree_results[i]
  phyloscanner_dir <-  str_replace(iqtree_results[i], "iqtree", "phyloscanner_results")

  if(!file.exists(phyloscanner_dir)){
    input_files <- c(input_files, iqtree_results1)
  }

}

input_files <- unlist(input_files)
all_inputs <- list()

for(j in 1:length(input_files)){
  all_input <- list.files(input_files[j], pattern="*.treefile", full.names = TRUE)
  all_inputs <- c(all_inputs, all_input)

}

all_inputs <- unlist(all_inputs)


#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(all_inputs, size = 10000)


mkdir_outfiles <- str_split(all_inputs, pattern = "/")
mkdir_outfiles <- unlist(lapply(mkdir_outfiles, function(x) paste(x[1:10], collapse = "/")))
mkdir_outfiles <- paste(mkdir_outfiles, "phyloscanner_results", sep = "/")

mkdir_split_dirs <- split_dirs(mkdir_outfiles, size = 10000)


#save input list file

#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

#input_split_dirs <- split_dirs(input_files, size = 10000)


for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_files_phyloscanner", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  #create files with mkdir file list
  mkdir_filename <- paste(paste("mkdir_files_phyloscanner", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = mkdir_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)



  mknew_dir <- paste("$(head -n $PBS_ARRAY_INDEX", mkdir_filename, "| tail -1)", sep = " ")


  pbsfilename <- paste(paste("phyloscanner", i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_dirs[[i]]), sep = "")




  pbstext <- paste("#PBS -l walltime=1:00:00",
                   "#PBS -l select=1:ncpus=1:mem=2gb",
                   paste("#PBS -o part2.stdout", sep = ""),
                   paste("#PBS -e part2.stderr", sep = ""),
                   array_number,
                   sep = "\n")

  pbstext <- paste(pbstext,
                   "\n",
                   "## load in the R environment",
                   "module load anaconda3/personal",
                   "source activate phyloscanner2",

                   "\n",
                   "## copy required files to the temporary directory on the compute node",
                   "export JOB_NUM=$(echo ${PBS_JOBID} | cut -f 1 -d '.' | cut -f 1 -d '[')",
                   "export WORKDIR=\"${EPHEMERAL}/${JOB_NUM}.${PBS_ARRAY_INDEX}\"",
                   "mkdir -p $WORKDIR",
                   "cd $WORKDIR",

                   "\n",
                   "#copy R scripts",
                   "cp $HOME/Ethics-HIV/large_pop/Rscripts/deepsequencing/8.Run_phyloscanner_part2_v4.R .",

                   paste("cp $HOME/Ethics-HIV/large_pop/deep_sequencing/villageModel/W0.01/phyloscanner_part2/750migrants/", input_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/large_pop/deep_sequencing/villageModel/W0.01/phyloscanner_part2/750migrants/", mkdir_filename, " .", sep = ""),

                   sep = "\n")



  pbstext <- paste(pbstext,

                   "\n",
                   "#create dir to copy treefiles to run phyloscanner",
                   "mkdir -p iqtree",
                   "#cp aligned treefiles to run phyloscanner",
                   paste(paste("cp -a ", "$(head -n $PBS_ARRAY_INDEX ", input_filename, "| tail -1)", sep = ""), "iqtree", sep = " "),

                   "\n",
                   "## Run R scripts",
                   "Rscript 8.Run_phyloscanner_part2_v4.R",

                   "\n",
                   "## make dir in the submission directory - anything not copied back will be lost",
                   paste(paste("mkdir -p ", "$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),


                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a ", "$WORKDIR/iqtree/*.csv", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),
                   paste("cp -a ", "$WORKDIR/iqtree/*.pdf", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),

                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}
