# create pbs script to run iqtree analysis on fasta files for reads
# to analyse with phyloscanner.
library(DescTools)
library(stringr)
library(HIVepisimAnalysis)

#reps that have to be completed for 500migrants
line_number <- 30
total_number_phyloscanner_windows <- 49

shiver_results <- list.files(list.files(list.files("/rds/general/user/fferre15/ephemeral/deepseq/best_trajectories_500migrants", full.names=T), full.names=T), full.names=T, pattern="shiver_results")
shiver_results <- shiver_results[line_number:length(shiver_results)]

#I initially run the job for phyloscanner to do all the
#rep each line in shiver_results to the total number of windows to analyse with
#phyloscanner
input_files <- shiver_results[sort(rep(seq_len(length(shiver_results)), total_number_phyloscanner_windows))]

#input_files <- paste(iqtree_results, "*.treefile", sep = "/")


#separate input files into chunks
#an array job can have a maximum of 10,000 subjobs

input_split_dirs <- split_dirs(input_files, size = 10000)


mkdir_outfiles <- str_split(input_files, pattern = "/")
mkdir_outfiles <- unlist(lapply(mkdir_outfiles, function(x) paste(x[1:10], collapse = "/")))
mkdir_outfiles <- paste(mkdir_outfiles, "phyloscanner_ali_340", sep = "/")

mkdir_split_dirs <- split_dirs(mkdir_outfiles, size = 10000)

window_numbers <- rep(1:total_number_phyloscanner_windows, length(input_files))
window_numbers_dirs <- split_dirs(window_numbers, size = 10000)


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

  #save window number
  windows_filename <- paste(paste("windows_phyloscanner", i, sep ="_"),
                            "txt", sep = ".")
  write.table(x = window_numbers_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = windows_filename)



  mknew_dir <- paste("$(head -n $PBS_ARRAY_INDEX", mkdir_filename, "| tail -1)", sep = " ")


  pbsfilename <- paste(paste("phyloscanner", i, sep = "_"), ".pbs", sep = "")


  array_number <- paste("#PBS", " -J", " 1-", length(input_split_dirs[[i]]), sep = "")




  pbstext <- paste("#PBS -l walltime=24:00:00",
                   "#PBS -l select=1:ncpus=1:mem=10gb",
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
                   "cp $HOME/Ethics-HIV/large_pop/Rscripts/deepsequencing/6.Run_phyloscanner_read_alignments_v8_340window.R .",
                   "cp $HOME/Ethics-HIV/large_pop/input_files/HIV1_REF_2020_phyloscanner.fasta .",
                   "#copy phyloscanner",
                   "cp $HOME/Ethics-HIV/Programs/phyloscanner.tar.gz .",

                   "\n",
                   paste("cp $HOME/Ethics-HIV/large_pop/deep_sequencing/villageModel/W0.01/phyloscanner_part1/500migrants/", input_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/large_pop/deep_sequencing/villageModel/W0.01/phyloscanner_part1/500migrants/", mkdir_filename, " .", sep = ""),
                   paste("cp $HOME/Ethics-HIV/large_pop/deep_sequencing/villageModel/W0.01/phyloscanner_part1/500migrants/", windows_filename, " .", sep = ""),

                   "\n",
                   "#untar phyloscanner",
                   "tar -xzvf phyloscanner.tar.gz",
                   sep = "\n")



  pbstext <- paste(pbstext,


                   "\n",
                   "## Run R scripts",
                   "Rscript 6.Run_phyloscanner_read_alignments_v8_340window.R \"/rds/general/user/fferre15/ephemeral/deepseq/best_trajectories_500migrants\" $PBS_ARRAY_INDEX",


                   "\n",
                   "## copy files back to Ephemeral directory",
                   "## I cannot copy results back to $HOME because of storage issues",
                   paste("cp -a $WORKDIR/AlignedReads", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),
                   paste("cp -a $WORKDIR/Consensuses", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),
                   paste("cp -a $WORKDIR/DuplicationData", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),
                   paste("cp -a $WORKDIR/ReadNames", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),
                   paste("cp -a $WORKDIR/WindowCoordinateCorrespondence.csv", paste("$(head -n $PBS_ARRAY_INDEX ", mkdir_filename, "| tail -1)", sep = ""), sep = " "),

                   "\n",
                   "rm -rf $WORKDIR",


                   sep = "\n")


  write.table(pbstext, file = pbsfilename, quote = FALSE,
              row.names = FALSE, col.names = FALSE)


}
