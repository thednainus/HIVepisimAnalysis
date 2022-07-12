#start of script
start_time <- Sys.time()

# You have to download IQ-TREE to run this script
# Change to the correct path of IQ-TREE on your computer
Software <- "/Applications/iqtree-2.1.2-MacOSX/bin/iqtree2"
#Software <- "iqtree"
maxCPU <- 2


#list files
#fasta_file <- list.files(".", pattern = "*.fasta")
#get argument on location of fasta file to read data
fasta_files <- list.files("phyloscanner_results/AlignedReads/", pattern = "fasta", full.names = T)
#fasta_file <- list.files(pattern = ".fasta")

for(i in 1:length(fasta_files)){
  ali_name <- paste("-s", fasta_files[i], sep = " ")

  #parameters to run iqtree
  iqtree_param <- c(ali_name, "-m HKY", "-T AUTO", "-ntmax", maxCPU)
  #run iqtree
  system2(command = Software, args = iqtree_param)
}




#start of script
end_time <- Sys.time()

iqtree_analyses <- data.frame(start = start_time, end = end_time)
saveRDS(iqtree_analyses, "iqtree_time.RDS")
