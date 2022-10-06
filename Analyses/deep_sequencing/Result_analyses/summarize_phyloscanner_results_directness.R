#summarize results from phyloscanner
#based on the whether two pairs are considered linked (independent of directness)
library(ape)
library(stringr)
library(lubridate)
library(dplyr)
library(HIVepisimAnalysis)

#beginning of simulations
init_sim_date <- ymd("1980-01-01")

#get all tips that was analysed
#I used a threshold for W of 0.55
#threshold <- 0.80

phyloscanner_results <- data.frame()

#common_dir <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/deepseq/W0.01"
common_dir <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/deepseq"

all_dirs <- list.files(list.files(list.files(common_dir, full.names = TRUE, pattern = "best_trajectories_*"),
                                  full.names = TRUE),
                       full.names = TRUE)

#all_dirs <- all_dirs[1:30]
#all_dirs <- all_dirs[1]

phyloscanner_results_dir <- "phyloscanner_results"

all_phyloscanner_results <- tibble()

reps_not_analysed <- c(31:100)

for(i in 1:length(all_dirs)){



  phyloscanner_dir <- paste(all_dirs[i], phyloscanner_results_dir, sep = "/")
  #print(i)
  #print(all_dirs[i])

  params <- str_split(phyloscanner_dir, pattern = "/")

  if(params[[1]][10] == "W0.01" ){
    mig <- str_split(params[[1]][11], pattern = "_")[[1]][3]
    param <- str_split(params[[1]][12], pattern = "_")[[1]][2]
    rep <- str_split(params[[1]][13], pattern = "_")[[1]][2]
  } else {
    mig <- str_split(params[[1]][10], pattern = "_")[[1]][3]
    param <- str_split(params[[1]][11], pattern = "_")[[1]][2]
    rep <- str_split(params[[1]][12], pattern = "_")[[1]][2]
  }

  if(!rep %in% reps_not_analysed){
      print(rep)
      filename <- "processing_network_results.tar.gz"


    dirs <- str_split(all_dirs[i], pattern = "/")


    if(dirs[[1]][10] == "W0.01" ){
      dir_location <- paste(str_split(all_dirs[i], pattern = "/")[[1]][c(1:9,11:13)], collapse = "/")
      filename <- paste(dir_location, filename, sep = "/")
      print(filename)

      #load transmission matrix
      untar(filename)
    }else{
      #load transmission matrix
      untar(paste(all_dirs[i], filename, sep = "/"))
    }

    load("output_deepseq/vts/merged_trees_sampling_migrant_years_1_simple_0.9.RData")

    #rep 1 , params 1067 (example for 1 replicate only)
    results <- read.csv(paste(phyloscanner_dir, "results_hostRelationshipSummary.csv", sep = "/"))

    #transmission phyloscanner
    trans_phy <- results %>%
      group_by(host.1, host.2) %>%
      group_modify(~ summarize_trans(.x))

    trans_phy["byGroup"] <- 1:nrow(trans_phy)


    #object "phylo_pairs" contain the information of whether two IDs represent
    #a transmission pair independent of the correct direction of transmission
    phylo_pairs <- trans_phy %>%
      group_by(byGroup) %>%
      group_modify(~check_true_pair(.x, tm))

    # check whether the pair host.1 and host.2 are a transmission pair
    # including the correct direction of transmission

    true_false_pairs_phylo <- phylo_pairs %>%
      group_by(byGroup) %>%
      group_modify(~check_true_transmissions(.x, tm))


    #here I will get all pairs that phyloscanner analysed and
    #add the sampled time for donor and recipient.
    pairs_with_dates <- true_false_pairs_phylo %>%
      group_by(byGroup) %>%
      group_modify(~add_sampled_times(.x, st_ids_region))


    #check for potential swaps of donor and recipients
    #that are the cases phyloscanner did directness == "yes"
    #real_pair == "yes" and trans =="false"
    pairs_with_dates["swaps"] <- "no"
    pairs_with_dates[pairs_with_dates$directness == "yes" &
                       pairs_with_dates$real_pair == "yes" &
                       pairs_with_dates$trans == "false",]$swaps <- "yes"

    all_trans_type <- pairs_with_dates

    all_trans_type["mig"] <- mig
    all_trans_type["param"] <- param
    all_trans_type["rep"] <- rep

    all_trans_type["total_pairs_analysed"] <- nrow(pairs_with_dates)
    all_trans_type["true_pairs_directness"] <- nrow(pairs_with_dates[pairs_with_dates$trans == "true",])
    all_trans_type["true_pairs_direction"] <- nrow(pairs_with_dates[pairs_with_dates$real_pair == "yes",])

    all_trans_type <- add_observed_values_directness(all_trans_type)

    all_phyloscanner_results <- rbind(all_phyloscanner_results, all_trans_type)

  }

  unlink("W_labels*")
  unlink("W_stats.csv")

}

#I will use the script below to create a function to sumarize the data
#these are results for direction only
saveRDS(all_phyloscanner_results, "all_phyloscanner_results_test_W80_directness.RDS")
saveRDS(all_phyloscanner_results, "all_phyloscanner_results_test_W0.01_directness.RDS")

#filtering by W > 0.80
all_phyloscanner_results_W80 <- readRDS("all_phyloscanner_results_test_W80_directness.RDS")
all_phyloscanner_results_W80["group"] <- paste(all_phyloscanner_results_W80$param,
                                               all_phyloscanner_results_W80$mig,
                                           sep = "_")
all_phyloscanner_results_W80["reps_groups"] <- paste(all_phyloscanner_results_W80$mig,
                                                 all_phyloscanner_results_W80$param,
                                                 all_phyloscanner_results_W80$rep,
                                                 sep = "_")

results_all_phylo_W80 <- all_phyloscanner_results_W80 %>%
  group_by(group) %>%
  mutate(precision = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FP")),
         sensitivity = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FN")),
         specificity = sum(observed == "TN")/(sum(observed == "TN") + sum(observed == "FP")),
         total_TP = sum(observed == "TP"),
         total_FP = sum(observed == "FP"),
         total_TN = sum(observed == "TN"),
         total_FN = sum(observed == "FN")) %>%
  select(group, precision, sensitivity, specificity,
         total_TP, total_FP, total_TN, total_FN) %>%
  distinct()

total_pairs_W80 <- all_phyloscanner_results_W80 %>%
  group_by(reps_groups, group) %>%
  distinct(total_pairs_analysed) %>%
  group_by(group) %>%
  mutate(total_pairs = sum(total_pairs_analysed)) %>%
  select(group, total_pairs) %>%
  distinct()



#filtering by W > 0.01
all_phyloscanner_results_W0.01 <- readRDS("all_phyloscanner_results_test_W0.01_directness.RDS")
all_phyloscanner_results_W0.01["group"] <- paste(all_phyloscanner_results_W0.01$param,
                                                 all_phyloscanner_results_W0.01$mig,
                                               sep = "_")
all_phyloscanner_results_W0.01["reps_groups"] <- paste(all_phyloscanner_results_W0.01$mig,
                                                       all_phyloscanner_results_W0.01$param,
                                                       all_phyloscanner_results_W0.01$rep,
                                                     sep = "_")

results_all_phylo_W0.01 <- all_phyloscanner_results_W0.01 %>%
  group_by(group) %>%
  mutate(precision = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FP")),
         sensitivity = sum(observed == "TP")/(sum(observed == "TP") + sum(observed == "FN")),
         specificity = sum(observed == "TN")/(sum(observed == "TN") + sum(observed == "FP")),
         total_TP = sum(observed == "TP"),
         total_FP = sum(observed == "FP"),
         total_TN = sum(observed == "TN"),
         total_FN = sum(observed == "FN")) %>%
  select(group, precision, sensitivity, specificity,
         total_TP, total_FP, total_TN, total_FN) %>%
  distinct()



total_pairs_W0.01 <- all_phyloscanner_results_W0.01 %>%
  group_by(reps_groups, group) %>%
  distinct(total_pairs_analysed) %>%
  group_by(group) %>%
  mutate(total_pairs = sum(total_pairs_analysed)) %>%
  select(group, total_pairs) %>%
  distinct()









