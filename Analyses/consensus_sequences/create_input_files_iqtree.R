#create input files for the trees that are missing
library(stringr)
library(HIVepisimAnalysis)

#get arguments
arguments <- commandArgs(trailingOnly = TRUE)



cluster <- "/rds/general/user/fferre15/ephemeral"
cluster_dir <- as.character(arguments[1])
#cluster_dir <- "best_trajectories_50migrants/params_2348"
dir_name <- paste(cluster, cluster_dir, sep = "/")

dirs <- dir(dir(dir_name, full.names = TRUE), full.names = TRUE, pattern = "perc*")

dirs1000bp <- paste(dirs, "output/vts/alignments/iqtree_results_1000bp", sep = "/")
dirs10000bp <- paste(dirs, "output/vts/alignments/iqtree_results_10000bp", sep = "/")


file_lists_1000bp <- list()
file_lists_10000bp <- list()

ml_file <- NULL
ml_file_10000bp <- NULL

for(i in 1:length(dirs1000bp)){

  #check whether directory exists
  if(file.exists(dirs1000bp[i]) == TRUE){

    #then check whether ML tree exists

    ml_file <- list.files(dirs1000bp[i], full.names = TRUE, pattern = "*.mldist")

  }

  if(length(ml_file) == 0 | file.exists(dirs1000bp[i]) == FALSE){

      #the list alignment to be used as input file to build tree with iqtree

      ali_name <- str_split(dirs1000bp[i], "/")
      ali_name <- unlist(lapply(ali_name, function(x) paste(x[1:13], collapse = "/")))
      fasta_file <- list.files(ali_name, full.names = TRUE, pattern = "*ali_1000bp")

      file_lists_1000bp <- c(file_lists_1000bp, fasta_file)


  }

}




for(j in 1:length(dirs10000bp)){

  #check whether directory exists
  if(file.exists(dirs10000bp[j]) == TRUE){

    #then check whether ML tree exists

    ml_file_10000bp <- list.files(dirs10000bp[j], full.names = TRUE, pattern = "*.mldist")

  }

  if(length(ml_file_10000bp) == 0 | file.exists(dirs10000bp[j]) == FALSE){

    #the list alignment to be used as input file to build tree with iqtree

    ali_name_10000bp <- str_split(dirs10000bp[j], "/")
    ali_name_10000bp <- unlist(lapply(ali_name_10000bp, function(x) paste(x[1:13], collapse = "/")))
    fasta_file_10000bp <- list.files(ali_name_10000bp, full.names = TRUE, pattern = "*ali_10000bp")

    file_lists_10000bp <- c(file_lists_10000bp, fasta_file_10000bp)


  }

}

#split dirs to have 10000 jobs per file
n = 10000
ali_1000bp <- unlist(file_lists_1000bp)

if(length(ali_1000bp) !=0){


  write.table(ali_1000bp, file = "ali_1000bp.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  input_split_dirs_1000bp <- split_dirs(ali_1000bp, size = n)

  for (i in 1:length(input_split_dirs_1000bp)){

    #create files with input file list
    input_filename <- paste(paste("input_file_list1000bp", i, sep ="_"),
                            "txt", sep = ".")
    write.table(x = input_split_dirs_1000bp[[i]], quote = FALSE,
                row.names = FALSE, col.names = FALSE,
                file = input_filename)


  }

}

ali_10000bp <- unlist(file_lists_10000bp)

if(length(ali_10000bp) !=0){

  write.table(ali_10000bp, file = "ali_10000bp.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  input_split_dirs_10000bp <- split_dirs(ali_10000bp, size = n)

  for (i in 1:length(input_split_dirs_10000bp)){

    #create files with input file list
    input_filename <- paste(paste("input_file_list10000bp", i, sep ="_"),
                            "txt", sep = ".")
    write.table(x = input_split_dirs_10000bp[[i]], quote = FALSE,
                row.names = FALSE, col.names = FALSE,
                file = input_filename)

  }
}
