#generate input dir files
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

ephemeral_filenames <- lapply(mkdir_names, function(x)
  paste("/rds/general/user/fferre15/ephemeral/results_small_pop", x, "results_sim.RDS", sep = "/"))

where2save <- lapply(mkdir_names, function(x)
  paste("/rds/general/user/fferre15/ephemeral/results_small_pop", x, sep = "/"))

#save filenames
for (i in 1:length(ephemeral_filenames)){

  #create files with input file list
  input_filename <- paste(paste("input_files_log_weights", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = ephemeral_filenames[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  where2save_filename <- paste(paste("where2save_files", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = where2save[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = where2save_filename)


  pbsfilename <- paste(paste("logweights_epidemics", i, sep = "_"), ".pbs", sep = "")
  array_number <- paste("#PBS", " -J", " 1-", length(ephemeral_filenames[[i]]), sep = "")


  pbstext <- paste("#PBS -l walltime=72:00:00",
                   "#PBS -l select=1:ncpus=1:mem=2gb",
                   "#PBS -o sim1_logweights.stdout",
                   "#PBS -e sim1_logweights.stderr",
                   array_number,
                   sep = "\n")

  pbstext <- paste(pbstext,
                   "\n",
                   "## load in the R environment",
                   "module load anaconda3/personal",
                   "source activate HIVepimodel",

                   "\n",
                   "## copy required files to the temporary directory on the compute node",
                   "export JOB_NUM=$(echo ${PBS_JOBID} | cut -f 1 -d '.' | cut -f 1 -d '[')",
                   "export WORKDIR=\"${EPHEMERAL}/${JOB_NUM}.${PBS_ARRAY_INDEX}\"",
                   "mkdir -p $WORKDIR",
                   "cd $WORKDIR",

                   "\n",
                   "#copy R scripts",

                   "cp $HOME/Ethics-HIV/small_pop/Rscripts/summarize_logs.R .",
                   paste("cp $HOME/Ethics-HIV/small_pop/MSM_population/input_files/", input_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/small_pop/MSM_population/input_files/", where2save_filename, " .", sep = ""),


                   "\n",
                   "## Run R using command line tags",
                   paste("Rscript summarize_logs.R $(head -n $PBS_ARRAY_INDEX ", input_filename, " | tail -1)", sep = ""),

                   "\n",
                   "## copy files back to the submission directory - anything not copied back will be lost",
                   paste("cp -a $WORKDIR/log_importance_wghts.RDS ", "$(head -n $PBS_ARRAY_INDEX ", where2save_filename, " | tail -1)", sep = ""),
                   sep = "\n")

  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)

}

