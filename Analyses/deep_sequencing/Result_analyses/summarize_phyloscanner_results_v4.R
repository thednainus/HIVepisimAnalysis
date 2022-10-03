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


    #get all pairs in which W >= 0.80
    #W_80 <- W1[W1$infectorProbability >= threshold,]

    #tips <- unique(c(W_80$donor, W_80$recip))
    #tips_tm <- unlist(lapply(tips, function(x) str_split(x, pattern = "_")[[1]][1]))
    #tips_id <- paste("ID", tips_tm, sep = "_")


    #trans_by_W80: all the real transmission that happened in the tm with correct direction
    #of transmission and W>= threshold.
    #trans_by_W80 <- real_trans[real_trans$infectorProbability >= threshold,]
    #names(trans_by_W80)[1:2] <- c("host.1", "host.2")
    #trans_by_W80$host.1 <- paste("ID", trans_by_W80$host.1, sep = "_")
    #trans_by_W80$host.2 <- paste("ID", trans_by_W80$host.2, sep = "_")

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

    #swap_pairs_phylo <- check_swap_transmissions(pairs_with_dates, tm)


    #tibble of true and swap pairs
    #true_swap <- rbind(true_pairs_phylo, swap_pairs_phylo)


    #remove from tibble trans_phy (transmissions by phyloscanner)
    #and check whether pairs represents real linked pairs of transmissions
    #pairs that still need to be analysed
    #to_analyse <- anti_join(phylo_pairs, true_swap,
    #                        by = c("host.1", "host.2"))

    #check for true pairs that has no ancestry
    #true_pairs_butNoAncestry <-    semi_join(to_analyse,
    #                                         trans_by_W80, by = c("host.1", "host.2"))

    #if(nrow(true_pairs_butNoAncestry) > 0){
    #  print("Check the below")
    #  print(i)
    #  print(all_dirs[i])

    #}


    # if(nrow(to_analyse) > 0){
    #
    #   if(nrow(true_pairs_butNoAncestry) == 0){
    #
    #     #to_analyse["order"] <- 1:nrow(to_analyse)
    #     linked_trans <- to_analyse %>%
    #       group_by(byGroup) %>%
    #       group_modify(~check_linked_transmissions(.x, tm))
    #     linked_trans <- linked_trans[,2:11]
    #     linked_trans["st_donor_ID"] <- NA
    #     linked_trans["st_recip_ID"] <- NA
    #
    #     all_trans_type <- rbind(true_swap[,c(2:13)],
    #                             linked_trans)
    #
    #   } else {
    #
    #     all_trans_type <- true_swap
    #
    #   }
    #
    #   all_trans_type["mig"] <- mig
    #   all_trans_type["param"] <- param
    #   all_trans_type["rep"] <- rep
    #
    # }

    all_trans_type <- pairs_with_dates

    all_trans_type["mig"] <- mig
    all_trans_type["param"] <- param
    all_trans_type["rep"] <- rep

    all_trans_type["total_pairs_analysed"] <- nrow(pairs_with_dates)
    all_trans_type["true_pairs_directness"] <- nrow(pairs_with_dates[pairs_with_dates$trans == "true",])
    all_trans_type["true_pairs_direction"] <- nrow(pairs_with_dates[pairs_with_dates$real_pair == "yes",])

    all_trans_type <- add_observed_values_direction(all_trans_type)

    all_phyloscanner_results <- rbind(all_phyloscanner_results, all_trans_type)

  }

  unlink("W_labels*")
  unlink("W_stats.csv")

}

#I will use the script below to create a function to sumarize the data

saveRDS(all_phyloscanner_results, "all_phyloscanner_results_test_W80.RDS")

all_phyloscanner_results <- readRDS("all_phyloscanner_results_test3.RDS")
all_phyloscanner_results["group"] <- paste(all_phyloscanner_results$param,
                                           all_phyloscanner_results$mig,
                                           sep = "_")
all_phyloscanner_results["reps_groups"] <- paste(all_phyloscanner_results$mig,
                                                 all_phyloscanner_results$param,
                                                 all_phyloscanner_results$rep,
                                                 sep = "_")

results_all_phylo <- all_phyloscanner_results %>%
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




results_all_phylo2 <- all_phyloscanner_results %>%
  group_by(group) %>%
  group_modify(~ summarize_all_data(.x))



direction_only["group"] <- paste(direction_only$param,
                                 direction_only$mig,
                                           sep = "_")
direction_only["reps_groups"] <- paste(direction_only$mig,
                                       direction_only$param,
                                       direction_only$rep,
                                                 sep = "_")

results_all_phylo <- direction_only %>%
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


summarize_all_data(all_phyloscanner_results)


