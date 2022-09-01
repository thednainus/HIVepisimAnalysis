#the change I made in this script is to rerun only the pairs that are missing from
#my previous analysis using a W threshold of 0.80
#my new threshold will be 0.01

# Simulate sequence alignments
# Seq-Gen MUST be installed in your computer
# http://tree.bio.ed.ac.uk/software/seqgen

# After computing W and deciding on a threshold, in this script
# I only keep in the tree, the tips in which pairs showed a W >= threshold.

#start of script
start_time <- Sys.time()


library(stringr)
library(ape)
library(seqinr)

#thrshold value for W (infector probability)
#for the time being this threshold will be 0.6
threshold1 <- 0.8
threshold2 <- 0.01

# Location for Seq-Gen. It should be changed to the correct location on your computer.
Software <- "/Applications/Seq-Gen-1.3.4/source/seq-gen"
#Software <- "Seq-Gen-1.3.4/source/seq-gen"
#reference genomic sequence to give as input to Seg-Gen
ref_genome <- "~/Desktop/Imperial/newHIVproject-01Aug2020/reference_subtypeB/genome/reference_subtypeB_completeGenome.phy"
#ref_genome <- "reference_subtypeB_completeGenome.phy"
# file containing the number of trees in which to simulate alignments
number_trees <- "~/Desktop/Imperial/newHIVproject-01Aug2020/reference_subtypeB/genome/tree_number.txt"
#number_trees <- "tree_number.txt"
# parameter for Seq-Gen
# seq-gen -mHKY -t0.5 -fe -n1 -k 1 < filename > filename+'.txt'
# -k parameter is to simulate sequences based on ancestral sequences.
# I will use HXB2 sequence to simulate HIV genomic sequences
# values for generating sequence alignment for HIV pol as described in van der Kyl and Berkhout, 2012
# values of -s0.0028 based on paper by Patino-Galindo and Gonzalez-Candelas 2017
# values of -t8.75 based on Chen et al 2004
# Updated parameters values based on sequences for San Diego + background sequences

#untar file
untar("processing_network_results.tar.gz")

#load file with infector probability information and phylogenetic trees
w_and_trees_filename <- list.files(path = "output_deepseq/vts",
                                   pattern ="*.RData", full.names = TRUE)
w_and_trees <- load(w_and_trees_filename)

#get tips that show W >= threshold
W_threshold1 <- W1[W1$infectorProbability >= threshold1,]
W_threshold2 <- W1[W1$infectorProbability >= threshold2,]

tipsToAnalyse1 <- unique(c(W_threshold1$donor, W_threshold1$recip))
tipsToAnalyse2 <- unique(c(W_threshold2$donor, W_threshold2$recip))

tipsToAnalyse <- tipsToAnalyse2[! tipsToAnalyse2 %in% tipsToAnalyse1]

# get all VirusTreeSimulator in which branch lengths have been converted into years
# and tip names are in the format of ID_migrant
vts_trees_years <- list.files(path = "output_deepseq/vts/merged_trees",
                              pattern ="*sampling.tre", full.names = TRUE)

tree_name <- vts_trees_years

filename_prefix <- str_split(tree_name, pattern = "\\.")[[1]][1]
filename_prefix <-str_split(filename_prefix, pattern = "/")[[1]][4]

# merge files with reference genome sequence and tree sequence for output
# to Seq_gen
if (!dir.exists("output_deepseq/vts/merged_trees/input_SeqGen")) {
  dir.create("output_deepseq/vts/merged_trees/input_SeqGen")
}

input_seqgen <- paste("output_deepseq/vts/merged_trees/input_SeqGen",
                      "tmp_inputSeqGen.txt",
                      sep = "/")


#read phylogenetic tree
vts_tree <- read.tree(tree_name)

#remove tips from the tree that will not be analysed with phyloscanner
all_IDs <- grepl(pattern = paste(paste("ID", tipsToAnalyse, sep="_"), collapse="|"),
                 x = vts_tree$tip.label)
all_IDs <- setNames(all_IDs, vts_tree$tip.label)

toDrop <- all_IDs[all_IDs == FALSE]
pairs_tree <- drop.tip(vts_tree, names(toDrop))

if(length(pairs_tree$tip.label) != (length(tipsToAnalyse)*10)){

  stop("Something is wrong with the grepl function as it is most
       likely mathing the incorrect number of tips")
}

#write tree to file to be used to simulate genetic sequence alignment
write.tree(phy = pairs_tree, file = "pairs_tree.tre")

#concatenate data for seqgen
system2(command = "cat",
        args = c(ref_genome, number_trees, "pairs_tree.tre"),
        stdout = input_seqgen)


#Create directory for sequence alignment
if (!dir.exists("output_deepseq/vts/merged_trees/alignments")) {
  dir.create("output_deepseq/vts/merged_trees/alignments")
}



# simulate sequence alignment for complete HIV genome
seq_filename_genome <-  paste("output_deepseq/vts/merged_trees/alignments/",
                              filename_prefix,
                              "_ali_genome.fasta",
                              sep = "")



# simulate sequence alignment using Seq-Gen
#Simulate complete HIV genomes
args1 <-  c("-mHKY",
            "-t8.75",
            "-f0.389,0.165,0.228,0.218",
            "-s0.0028",
            "-n1",
            "-of",
            "-k 1",
            input_seqgen)

system2(command = Software,
        args = args1,
        stdout = seq_filename_genome)

# read alignment and split each combination of sequences per individual
# to simulate Illumina reads
genome_seqs <- read.FASTA(seq_filename_genome)

#get ID names
ID_names <-  unlist(lapply(names(genome_seqs),
                           function(x)str_split(string = x, pattern = "_")[[1]][2]))

#create dataframe with sequence name in phylogenetic tree and ID names
seq_ID_names <-  data.frame(seq_name = names(genome_seqs), ID_names = ID_names)
seq_ID_names$ID_names <- as.factor(seq_ID_names$ID_names)
seq_ID_split <- split(seq_ID_names, f = seq_ID_names$ID_names)

#match id_names to sequence name in DNAbin
index <- lapply(seq_ID_split, function(x) match(unname(unlist(x[[1]])),
                                                names(genome_seqs)))
seq_by_ID <- lapply(index, function(x)  genome_seqs[x])

#save separate alignment in fasta format

if (!dir.exists("output_deepseq/vts/merged_trees/alignments/by_ID")) {
  dir.create("output_deepseq/vts/merged_trees/alignments/by_ID")
}

dir_name_by_tree <- paste("output_deepseq/vts/merged_trees/alignments/by_ID",
                          filename_prefix,
                          sep = "/")

if (!dir.exists(dir_name_by_tree)) {
  dir.create(dir_name_by_tree)
}


lapply(seq_by_ID,
       function(x) write.FASTA(x,file = paste( dir_name_by_tree, "/",
                                               paste("ID_", str_split(names(x)[1],
                                                                      pattern = "_")[[1]][2],
                                                     ".fasta", sep = ""), sep = "")))

#read fasta files by ID and split it by amplicon size
#amplicon size will follow Gal et all. 2012
by_ID_ali <- list.files(dir_name_by_tree, pattern = "fasta",full.names = TRUE)

ampl_dirname <- paste(paste(str_split(by_ID_ali[1], "/")[[1]][1:6], collapse = "/"),
                      "amplicons", sep = "/")

if (!dir.exists(ampl_dirname)) {
  dir.create(ampl_dirname)
}

for (ali in by_ID_ali){
  #read alignment
  IDali <- read.fasta(ali)
  ampl1 <- lapply(IDali, function(x) x[480:2407])
  names(ampl1) <- paste(names(ampl1), "ampl1", sep = "_")
  ampl2 <- lapply(IDali, function(x) x[1487:5058])
  names(ampl2) <- paste(names(ampl2), "ampl2", sep = "_")
  ampl3 <- lapply(IDali, function(x) x[4783:7848])
  names(ampl3) <- paste(names(ampl3), "ampl3", sep = "_")
  ampl4 <- lapply(IDali, function(x) x[5968:9518])
  names(ampl4) <- paste(names(ampl4), "ampl4", sep = "_")

  all_ampl <- c(ampl1, ampl2, ampl3, ampl4)

  filename_ali <- str_split(ali, pattern = "/")[[1]][7]
  write.fasta(sequences = all_ampl, names = names(all_ampl),
              file.out = paste(ampl_dirname, filename_ali, sep = "/"))

}






#end of script
end_time <- Sys.time()
print("Sequence alignment simulation took:")
end_time - start_time

seq_alignment_time <- data.frame(start = start_time, end = end_time)
saveRDS(seq_alignment_time, "seq_alignment_time_deepseq.RDS")
