#delete Illumina_reads directory in the cluster if shiver analyses
#have been carried out
#This is to relase storage space as the Illumina reads files occupy a
#large amount of data

library(stringr)
library(HIVepisimAnalysis)


cluster <- "/ocean/projects/bio210082p/nascimen/data_deepsequencing"
dir_name <- paste(cluster, "best_trajectories_50migrants/params_1067", sep = "/")
dirs=dir(dir(dir_name, full.names = TRUE), full.names = TRUE, pattern = "perc*")

#count number of Illumina files
dirs_illumina <- paste(dirs, "Illumina_reads/results_sampling", sep = "/")

#list of Illumina dirs that do not have all equivalent shiver files
dir_list_to_do <- list()
ids_list_to_do <- list()

#do a loop and compare number of files with Illumina reads and shiver files

for(i in 1:length(dirs_illumina)){

  files_illumina <- dir(dirs_illumina[i], full.names = TRUE)
  shiver_files <- paste(dirs[i], "shiver_results", sep = "/")
  files_shiver <- dir(shiver_files, full.names = TRUE)

  if(((length(files_shiver)/7) == length(files_illumina)) & (length(files_illumina) != 0)){

    #rm Illumina directory

    command_rm <- paste("rm -rf", dirs_illumina[i])
    system(command_rm)



  }

  if(((length(files_shiver)/7) != length(files_illumina)) & (length(files_illumina) != 0)){

    dir_list_to_do = c(dir_list_to_do, dirs_illumina[i])

    illumina_ID_names <- dir(dirs_illumina[i], full.names = FALSE)
    files_shiver2 <- dir(shiver_files, full.names = FALSE, pattern = "*remap.bam.bai")

    #split names
    shiver_ID_names <- str_split(files_shiver2, "_")
    shiver_ID_names <- unlist(lapply(shiver_ID_names, function(x) paste(x[1:2], collapse = "_")))

    #get index of the ID that is still missing in shiver files
    index <- which(! illumina_ID_names %in% shiver_ID_names)

    #to do Illumina to shiver
    ids_list_to_do <- c(ids_list_to_do, files_illumina[index])


  }

}

to_do_list <- unlist(dir_list_to_do)
write.table(to_do_list, file = "toDo_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

illumina_to_do_list <- unlist(ids_list_to_do)
write.table(illumina_to_do_list, file = "IDs_toDo_list.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

#split dirs to have 1000 jobs per file
n = 1000

input_files <- paste(illumina_to_do_list, "/.", sep = "")
input_split_dirs <- split_dirs(input_files, size = n)


output_files <- str_split(input_files, pattern = "/")
output_files <- unlist(lapply(output_files, function(x) paste(x[1:10], collapse = "/")))
output_split_dirs <- split_dirs(output_files, size = n)

#mknew dir to copy shiver files
mkdir_outfiles <- lapply(output_split_dirs, function(x) paste(x, "shiver_results", sep = "/"))

for (i in 1:length(input_split_dirs)){

  #create files with input file list
  input_filename <- paste(paste("input_file_list", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = input_split_dirs[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = input_filename)

  mkdir_filename <- paste(paste("mkdir_file_list", i, sep ="_"),
                          "txt", sep = ".")
  write.table(x = mkdir_outfiles[[i]], quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = mkdir_filename)

}

