# Estimate phylogenetic tree using IQ-TREE
# to use with phyloscanner
# input file will be the read alignmnets created using script
# 6.Run_phyloscanner_read_alignments_v2.R

#start of script
start_time <- Sys.time()


# You have to download IQ-TREE to run this script
# Change to the correct path of IQ-TREE on your computer
Software <- "iqtree"
maxCPU <- 2


#list files
#fasta_file <- list.files(".", pattern = "*.fasta")
#get argument on location of fasta file to read data
fasta_file <- list.files(pattern = ".fasta")


ali_name <- paste("-s", fasta_file, sep = " ")

#parameters to run iqtree
iqtree_param <- c(ali_name, "-m HKY", "-T AUTO", "-ntmax", maxCPU)
#run iqtree
system2(command = Software, args = iqtree_param)


#start of script
end_time <- Sys.time()

iqtree_analyses <- data.frame(start = start_time, end = end_time)
saveRDS(iqtree_analyses, "iqtree_time.RDS")
