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

#distance threshold
distance <- "--distanceThreshold 0.05"

#normRefFileName: using
#update this later using the 2019 reference genome sequences I used with shiver
#check manual on how to create this normalization file
norm <- "--normRefFileName /Users/user/Desktop/Imperial/newHIVproject-01Aug2020/Phyloscanner/phyloscanner/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv"

# downsamplig option
#downsampling <- "--maxReadsPerHost 2"


# list directories
output_dirs <- dir(path = "output_deepseq/vts/merged_trees/Illumina_reads",
                   full.names = TRUE)
output_dirs <- output_dirs[1]

for(i in 1:length(output_dirs)){

  #directory phyloscanner used to save results of make_trees
  pyloscanner_fullPath <- paste(SplitPath(output_dirs[i])$normpath,
                                    "phyloscanner", sep = "/")

  #location of raxml trees
  raxml_trees <- paste(pyloscanner_fullPath,
                           "RAxMLfiles",
                           "RAxML_bestTree.InWindow_", sep = "/")


  parameters <- paste(raxml_trees, output_label, splitsRule, outgroup,
                      multifurcation, distance, norm, sep = " ")


  command_analyse_trees <- paste("cd", pyloscanner_fullPath, "&&", analyse_trees, sep = " ")
  analyseTrees_and_args <- paste(command_analyse_trees, parameters, sep = " ")


  system(analyseTrees_and_args)





}

#end of script
end_time <- Sys.time()
print("Simulation took:")
end_time - start_time
