#script to get the replicates that still needs to be analysed

#create input files for the trees that are missing
library(stringr)
library(HIVepisimAnalysis)

#get arguments
#arguments <- commandArgs(trailingOnly = TRUE)
arguments <- "best_trajectories_750migrants"



#cluster <- "/rds/general/user/fferre15/ephemeral"
cluster <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper"
#cluster_dir <- as.character(arguments[1])
cluster_dir <- arguments
#cluster_dir <- "best_trajectories_50migrants/params_2348"
dir_name <- paste(cluster, cluster_dir, sep = "/")

#sampler1
#dirs <- dir(dir(dir(dir_name, full.names = TRUE), full.names = TRUE),
#            full.names = TRUE, pattern = "perc*")

#sampler2
dirs <- dir(dir(dir(dir(dir_name, full.names = TRUE), full.names = TRUE),
            full.names = TRUE, pattern = "sampler2*"), full.names = TRUE, pattern = "perc*")

alignment_dir <- paste(dirs, "output/vts/alignments/W_estimated", sep = "/")


file_lists_1000bp <- data.frame()
file_lists_10000bp <- data.frame()

Wstats_file_1000bp <- NULL
Wstats_file_10000bp <- NULL

for(i in 1:length(alignment_dir)){

  #check whether directory exists
  if(file.exists(alignment_dir[i]) == TRUE){

    #then check whether ML tree exists

    Wstats_1000bp <- list.files(alignment_dir[i], full.names = TRUE, pattern = "1000bp_W_stats*")
    Wstats_10000bp <- list.files(alignment_dir[i], full.names = TRUE, pattern = "10000bp_W_stats*")

  }

  if(length(Wstats_1000bp) == 0 | file.exists(alignment_dir[i]) == FALSE){

    #the list alignment to be used as input file to build tree with iqtree

    ali_name <- str_split(alignment_dir[i], "/")
    #param <- str_split(ali_name[[1]][10], pattern = "_")[[1]][2]
    #rep <- str_split(ali_name[[1]][11], pattern = "_")[[1]][2]
    #perc <- str_split(ali_name[[1]][12], pattern = "_")[[1]][2]

    param <- str_split(ali_name[[1]][10], pattern = "_")[[1]][2]
    rep <- str_split(ali_name[[1]][11], pattern = "_")[[1]][2]
    perc <- str_split(ali_name[[1]][13], pattern = "_")[[1]][2]

    params <- data.frame(param = param, perc = perc, rep = rep)

    file_lists_1000bp <- rbind(file_lists_1000bp, params)


  }


  if(length(Wstats_10000bp) == 0 | file.exists(alignment_dir[i]) == FALSE){

    #the list alignment to be used as input file to build tree with iqtree

    ali_name <- str_split(alignment_dir[i], "/")
    #param <- str_split(ali_name[[1]][10], pattern = "_")[[1]][2]
    #rep <- str_split(ali_name[[1]][11], pattern = "_")[[1]][2]
    #perc <- str_split(ali_name[[1]][12], pattern = "_")[[1]][2]

    param <- str_split(ali_name[[1]][10], pattern = "_")[[1]][2]
    rep <- str_split(ali_name[[1]][11], pattern = "_")[[1]][2]
    perc <- str_split(ali_name[[1]][13], pattern = "_")[[1]][2]

    params <- data.frame(param = param, perc = perc, rep = rep)

    file_lists_10000bp <- rbind(file_lists_10000bp, params)


  }

}





