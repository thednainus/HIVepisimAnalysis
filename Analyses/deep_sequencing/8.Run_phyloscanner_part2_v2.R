# script to run phyloscanner: https://github.com/BDI-pathogens/phyloscanner
# this script run the second part of phyloscanner:
# it will analyze trees generated with script 7.Estimate_trees_read_alignments_v2.R

# phyloscanner dependencies
# https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh

#start of script
start_time <- Sys.time()

library(DescTools)
library(reticulate)

#location of phyloscanner make analyse_trees R program
phyloscanner <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner"
analyse_trees <- paste(phyloscanner, "phyloscanner_analyse_trees.R", sep = "/")

# phyloscanner analyse_trees parameters that does not depend on location
# of results of running run_phyloscanner_part1.R

#label to name output results
output_label <- "results_"

#splits rule
splitsRule <- "s,20"

# outgroup name
outgroup <- "--outgroupName C.BW.00.00BW07621.AF443088"

# multifurcation threshold as g = "guess"
multifurcation <- "--multifurcationThreshold g"

#allow multi transmission
multtrans <- "--allowMultiTrans"

confidenceCut <- "--confidence.cut 0.5"

#distance threshold
#summarize all paiwise relationships in csv summary file
# I did not used the value for distance threshold and phyloscanner then assumes
#it is Inf (shows all transmission pairs)
#distance <- "--distanceThreshold 0.05"

#normRefFileName: using
#update this later using the 2019 reference genome sequences I used with shiver
#check manual on how to create this normalization file
norm_file_location <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv"
norm <- paste("--normRefFileName", norm_file_location, sep = " ")

# downsamplig option
downsampling <- "--maxReadsPerHost 50"


# list directories
tree_files <- list.files(path = "/Users/user/Desktop/teste/iqtree/",
                   pattern = "treefile")

#directory phyloscanner used to save results of make_trees
iqtree_dir <- "/Users/user/Desktop/teste/iqtree"

#location of iqtree tree files
iqtree_trees <- paste(iqtree_dir, "*.treefile", sep = "/")

treefile_extension <- "--treeFileExtension .treefile"

parameters <- paste(iqtree_dir, treefile_extension, output_label, splitsRule, outgroup,
                    multifurcation, multtrans, norm, sep = " ")


command_analyse_trees <- paste("cd", iqtree_dir, "&&", analyse_trees, sep = " ")
analyseTrees_and_args <- paste(command_analyse_trees, parameters, sep = " ")


system(analyseTrees_and_args)



#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time
