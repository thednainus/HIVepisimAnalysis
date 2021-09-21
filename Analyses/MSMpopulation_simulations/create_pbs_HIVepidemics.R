library(stringr)
library(HIVepisimAnalysis)

file_numbers <- list()

for(i in 1:1000){

  file_numbers[[i]] <- rep(i, 50)

}

file_numbers <- unlist(file_numbers)

file_names <- split_dirs(file_numbers, size = 10000)

mkdir_file_names <- lapply(file_names, function(x) paste("small_msm_pop_sim",x, sep = "_"))

reps <- c(1:50)
reps_total <- rep(reps, 200)

mkdir_names <- lapply(mkdir_file_names, function(x) paste(x, paste("rep_", reps_total, sep = ""), sep = "/"))

#save filenames
for (i in 1:length(file_names)){

  #create files with input file list
  input_filename <- paste(paste("input_files_msmpopulation", "set_sim", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = file_names[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  #create files with mkdir filenames list
  mkdir_filename <- paste(paste("mkdir_files_msmpopulation", "set_sim", i, sep ="_"),
                          "txt", sep = ".")

  write.table(x = mkdir_names[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)


  pbsfilename <- paste(paste("msmpopulation_epidemics", "sim_set", i, sep = "_"), ".pbs", sep = "")
  array_number <- paste("#PBS", " -J", " 1-", length(file_names[[i]]), sep = "")


  pbstext <- paste("#PBS -l walltime=72:00:00",
                   "#PBS -l select=1:ncpus=1:mem=2gb",
                   "#PBS -o sim1_epidemics.stdout",
                   "#PBS -e sim1_epidemics.stderr",
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

                   "cp $HOME/Ethics-HIV/small_pop/Rscripts/Simulate_MSMpopulation_small_pop_lhs.R .",
                   "cp $HOME/Ethics-HIV/small_pop/MSM_population/input_files/input_files_msmpopulation_set_sim_1.txt .",
                   "cp $HOME/Ethics-HIV/small_pop/MSM_population/input_files/mkdir_files_msmpopulation_set_sim_1.txt .",

                   "\n",
                   "## Run R using command line tags",
                   "Rscript Simulate_MSMpopulation_small_pop_lhs.R $(head -n $PBS_ARRAY_INDEX input_files_msmpopulation_set_sim_1.txt| tail -1)",

                   "\n",
                   "## make dir in Ephemeral",
                   "mkdir -p $EPHEMERAL/$(head -n $PBS_ARRAY_INDEX mkdir_files_msmpopulation_set_sim_1.txt| tail -1)",

                   "\n",
                   "## copy files back to the submission directory - anything not copied back will be lost",
                   "cp -a $WORKDIR/. $EPHEMERAL/$(head -n $PBS_ARRAY_INDEX mkdir_files_msmpopulation_set_sim_1.txt| tail -1)",
                   sep = "\n")

  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)

}
