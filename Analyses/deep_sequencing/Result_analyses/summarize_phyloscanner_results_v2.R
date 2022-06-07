#summarize results from phyloscanner
library(ape)
library(stringr)
library(lubridate)
library(dplyr)

phyloscanner_results <- data.frame()

#get all tips that was analysed
#I used a threshold for W of 0.55
threshold <- 0.55

#get all pairs in which W >= 0.55
W_55 <- W1[W1$infectorProbability >= 0.55,]

tips <- unique(c(W_55$donor, W_55$recip))
tips_tm <- unlist(lapply(tips, function(x) str_split(x, pattern = "_")[[1]][1]))

#get all pairs involving any individuals in which W >= threshold
tm_subset_list <- lapply(tips_tm, function(x) tm[tm$sus == x | tm$inf == x ,])
tm_subset_df <- do.call(rbind, tm_subset)

#get sampling time of all individuals in tips
sample_times <- lapply(tips, function(x) st_ids_region[st_ids_region$tip_name == x,])
sample_times <- do.call(rbind, sample_times)

#estimate difference between sampling time and time of infection
get_difference <- function(df1, df2){
  #df1=dataframe of sampling times
  #df2=datafame of sampling time
}

lapply(sample_times$sampled_ID, function(x) tm_subset_df[tm_subset_df$sus == x,])

#beginning of simulations
init_sim_date <- ymd("1980-01-01")

common_dir <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants/params_1067/rep_1/iqtree"


#load transmission matrix
untar("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/best_trajectories_500migrants/params_1067/rep_1/processing_network_results.tar.gz")
load("output_deepseq/vts/merged_trees_sampling_migrant_years_1_simple_0.9.RData")
#real_trans_W <- readRDS("real_trans_W_tm.RDS")

trans_by_W <- real_trans
names(trans_by_W)[1:2] <- c("host.1", "host.2")
trans_by_W$host.1 <- paste("ID", trans_by_W$host.1, sep = "_")
trans_by_W$host.2 <- paste("ID", trans_by_W$host.2, sep = "_")

#rep 1 , params 1067 (example for 1 replicate only)
results <- read.csv(paste(common_dir, "results_rep_1__hostRelationshipSummary.csv", sep = "/"))
results_subset <- subset(results, ancestry != "noAncestry")

#correct donor and recipient
all_trans <- semi_join(results_subset, trans_by_W, by = c("host.1", "host.2"))
#semi_join(tm_region1, all_trans_5perc, by = c("host.1", "host.2"))
#convert fraction to decimal number
all_trans$fraction.math <- sapply(all_trans$fraction, function(x) eval(parse(text=x)))
#all_trans_subset <- subset(all_trans, ancestry != "noAncestry" & ancestry != "multiTrans")


#swap of donor and recipient
trans_by_W_trunc <- data.frame(host.1 = trans_by_W$host.2, host.2 = trans_by_W$host.1)
all_trans2 <- semi_join(results_subset, trans_by_W_trunc, by = c("host.1", "host.2"))
all_trans2$fraction.math <- sapply(all_trans2$fraction, function(x) eval(parse(text=x)))
all_trans2_subset <- subset(all_trans_5perc2, ancestry != "noAncestry" & ancestry != "multiTrans")


perc5 <- data.frame(total_pairs = nrow(sp5perc),
                    correct_pair = nrow(all_trans_5perc_subset),
                    swap_pairs = nrow(all_trans_5perc2_subset), perc = "5perc")

phyloscanner_results <- rbind(phyloscanner_results, perc5)



tm_region <- subset(tm, infOrigin == "region" & susOrigin == "region")
tm_region$year <- days2years(tm_region$at, init_date = init_sim_date)

#analysis of 5% -----
#get tm in which sus got infected before sampling time
#load tree for 5%

tree5perc <- read.tree("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/deep_sequencing_results_cluster/5perc/results_sampling.tre")
tree5perc$tip.label <- as.numeric(unlist(lapply(tree5perc$tip.label,
                                                function(x) str_split(x, "_")[[1]][2])))
#keep_only_rows that tips are in phylogenetic tree
rows_to_keep1 <- keep_row(df = tm_region, tree = tree5perc)
tm5perc <- tm_region[rows_to_keep1,]


#get transmission that happened before sampling time
tm5perc_bst <- apply(st_ids_region, 1, function(x) tm5perc[tm5perc$sus == x[1] & tm5perc$year <= x[2],])
tm5perc_bst <- do.call(rbind, tm5perc_bst)



tm_region1 <- tm5perc_bst
tm_region2 <- tm5perc_bst
tm_region1$host.1 <- paste("ID", tm5perc_bst$inf, sep = "_")
tm_region1$host.2 <- paste("ID", tm5perc_bst$sus, sep = "_")


tm_region2$host.1 <- paste("ID", tm5perc_bst$sus, sep = "_")
tm_region2$host.2 <- paste("ID", tm5perc_bst$inf, sep = "_")





#sampling proportion = 5% ----
sp5perc <- read.csv(paste(common_dir, "5perc/iqtree/treefile/phyloscanner_results/results_5perc__hostRelationshipSummary.csv", sep = "/"))
sp5perc_subset <- subset(sp5perc, ancestry != "noAncestry")

#correct donor and recipient
all_trans_5perc <- semi_join(sp5perc, tm_region1, by = c("host.1", "host.2"))
semi_join(tm_region1, all_trans_5perc, by = c("host.1", "host.2"))
#convert fraction to decimal number
all_trans_5perc$fraction.math <- sapply(all_trans_5perc$fraction, function(x) eval(parse(text=x)))
all_trans_5perc_subset <- subset(all_trans_5perc, ancestry != "noAncestry" & ancestry != "multiTrans")


#swap of donor and recipient
all_trans_5perc2 <- semi_join(sp5perc, tm_region2, by = c("host.1", "host.2"))
all_trans_5perc2$fraction.math <- sapply(all_trans_5perc2$fraction, function(x) eval(parse(text=x)))
all_trans_5perc2_subset <- subset(all_trans_5perc2, ancestry != "noAncestry" & ancestry != "multiTrans")


perc5 <- data.frame(total_pairs = nrow(sp5perc),
                    correct_pair = nrow(all_trans_5perc_subset),
                    swap_pairs = nrow(all_trans_5perc2_subset), perc = "5perc")

phyloscanner_results <- rbind(phyloscanner_results, perc5)


#sampling proportion = 30% ----


#get tm in which sus got infected before sampling time
#load tree for 30%

tree30perc <- read.tree("/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/deep_sequencing_results_cluster/30perc/results_sampling.tre")
tree30perc$tip.label <- as.numeric(unlist(lapply(tree30perc$tip.label,
                                                function(x) str_split(x, "_")[[1]][2])))
#keep_only_rows that tips are in phylogenetic tree
rows_to_keep2 <- keep_row(df = tm_region, tree = tree30perc)
tm30perc <- tm_region[rows_to_keep2,]


#get transmission that happened before sampling time
#load trasmission matrix
load(paste(common_dir, "30perc/merged_trees_sampling_migrant_years_1_simple_0.035.RData", sep = "/"))

tm30perc_bst <- apply(st_ids_region, 1, function(x) tm30perc[tm30perc$sus == x[1] & tm30perc$year <= x[2],])
tm30perc_bst <- do.call(rbind, tm30perc_bst)



tm_region3 <- tm30perc_bst
tm_region4 <- tm30perc_bst
tm_region3$host.1 <- paste("ID", tm30perc_bst$inf, sep = "_")
tm_region3$host.2 <- paste("ID", tm30perc_bst$sus, sep = "_")


tm_region4$host.1 <- paste("ID", tm30perc_bst$sus, sep = "_")
tm_region4$host.2 <- paste("ID", tm30perc_bst$inf, sep = "_")


sp30perc <- read.csv(paste(common_dir, "30perc/iqtree/results_30perc__hostRelationshipSummary.csv", sep = "/"))





#correct donor and recipient
all_trans_30perc <- semi_join(sp30perc, tm_region3, by = c("host.1", "host.2"))
all_trans_in_30perc <- semi_join(tm_region3, all_trans_30perc, by = c("host.1", "host.2"))
#convert fraction to decimal
all_trans_30perc$fraction.math <- sapply(all_trans_30perc$fraction, function(x) eval(parse(text=x)))
#subset data
all_trans_30perc_subset <- subset(all_trans_30perc, ancestry != "noAncestry" &
                                    ancestry != "multiTrans" & ancestry != "complex")


#swap of donor and recipient
all_trans_30perc2 <- semi_join(sp30perc, tm_region4, by = c("host.1", "host.2"))
all_trans_30perc2$fraction.math <- sapply(all_trans_30perc2$fraction, function(x) eval(parse(text=x)))
all_trans_30perc2_subset <- subset(all_trans_30perc2, ancestry != "noAncestry" &
                                     ancestry != "multiTrans" &
                                     ancestry != "complex")


perc30 <- data.frame(total_pairs = nrow(sp30perc),
                    correct_pair = nrow(all_trans_30perc_subset),
                    swap_pairs = nrow(all_trans_30perc2_subset), perc = "30perc")



#n_tm = number of true transmission pairs before sampling time
#total_pairs = total pairs analysed by phyloscanner
phyloscanner_results["n_tm"] <- nrow(tm5perc_bst)
perc30["n_tm"] <- nrow(tm30perc_bst)
phyloscanner_results <- rbind(phyloscanner_results, perc30)

saveRDS(phyloscanner_results, "phyloscanner_results.RDS")

