#summarize results from phyloscanner
library(ape)
library(stringr)
library(lubridate)
library(dplyr)
library(HIVepisimAnalysis)

#beginning of simulations
init_sim_date <- ymd("1980-01-01")

#get all tips that was analysed
#I used a threshold for W of 0.55
threshold <- 0.80

phyloscanner_results <- data.frame()

common_dir <- "/Users/user/Desktop/Imperial/newHIVproject-01Aug2020/R_projects/Results_paper/deepseq"

all_dirs <- list.files(list.files(list.files(common_dir, full.names = TRUE),
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
  mig <- str_split(params[[1]][10], pattern = "_")[[1]][3]
  param <- str_split(params[[1]][11], pattern = "_")[[1]][2]
  rep <- str_split(params[[1]][12], pattern = "_")[[1]][2]

  if(!rep %in% reps_not_analysed){
    print(rep)
    filename <- "processing_network_results.tar.gz"

    #load transmission matrix
    untar(paste(all_dirs[i], filename, sep = "/"))
    load("output_deepseq/vts/merged_trees_sampling_migrant_years_1_simple_0.9.RData")

    #get all pairs in which W >= 0.80
    W_80 <- W1[W1$infectorProbability >= threshold,]

    tips <- unique(c(W_80$donor, W_80$recip))
    tips_tm <- unlist(lapply(tips, function(x) str_split(x, pattern = "_")[[1]][1]))
    tips_id <- paste("ID", tips_tm, sep = "_")


    #trans_by_W80: all the real transmission that happened in the tm with correct direction
    #of transmission and W>= threshold.
    trans_by_W80 <- real_trans[real_trans$infectorProbability >= threshold,]
    names(trans_by_W80)[1:2] <- c("host.1", "host.2")
    trans_by_W80$host.1 <- paste("ID", trans_by_W80$host.1, sep = "_")
    trans_by_W80$host.2 <- paste("ID", trans_by_W80$host.2, sep = "_")

    #rep 1 , params 1067 (example for 1 replicate only)
    results <- read.csv(paste(phyloscanner_dir, "results_hostRelationshipSummary.csv", sep = "/"))

    #transmission phyloscanner
    trans_phy <- results %>%
      group_by(host.1, host.2) %>%
      group_modify(~ summarize_trans(.x))

    #check whether the direction of transmission in phyloscanner corresponds to
    #true transmission
    all_trans <- semi_join(trans_phy,
                           trans_by_W80, by = c("host.1", "host.2"))

    true_pairs_phylo <- check_true_transmissions(trans_phy, all_trans)

    #host.1 is the donor
    true_pairs_phylo["st_donor_ID"] <- unlist(lapply(true_pairs_phylo$host.1, function(x)
      st_ids_region[paste("ID", st_ids_region$sampled_ID, sep = "_") == x,]$sampled_time))

    #host.2 is the recipient
    true_pairs_phylo["st_recip_ID"] <- unlist(lapply(true_pairs_phylo$host.2, function(x)
      st_ids_region[paste("ID", st_ids_region$sampled_ID, sep ="_") == x,]$sampled_time))

    #time that transmission happened according to transmission matrix
    true_pairs_phylo["time_trans"] <- apply(true_pairs_phylo, 1, function(x)
      tm[tm$sus == str_split(x[[2]], pattern = "_")[[1]][2] &
           tm$inf == str_split(x[[1]], pattern = "_")[[1]][2],]$year)


    #swap of donor and recipient
    trans_by_W_trunc <- data.frame(host.1 = trans_by_W80$host.2,
                                   host.2 = trans_by_W80$host.1)
    all_trans2 <- semi_join(trans_phy, trans_by_W_trunc,
                            by = c("host.1", "host.2"))

    swap_pairs_phylo <- check_swap_transmissions(trans_phy, all_trans2)


    #tibble of true and swap pairs
    true_swap <- rbind(true_pairs_phylo, swap_pairs_phylo)


    #remove from tibble trans_phy (transmissions by phyloscanner)
    #and check whether pairs represents real linked pairs of transmissions
    #pairs that still need to be analysed
    to_analyse <- anti_join(trans_phy, true_swap,
                            by = c("host.1", "host.2"))

    #check for true pairs that has no ancestry
    true_pairs_butNoAncestry <-    semi_join(to_analyse,
                                             trans_by_W80, by = c("host.1", "host.2"))

    if(nrow(true_pairs_butNoAncestry) > 0){
      print("Checke the below")
      print(i)
      print(all_dirs[i])

    }


    if(nrow(to_analyse) > 0){

      if(nrow(true_pairs_butNoAncestry) == 0){

        to_analyse["order"] <- 1:nrow(to_analyse)
        linked_trans <- to_analyse %>%
          group_by(order) %>%
          group_modify(~check_linked_transmissions(.x, tm))
        linked_trans <- linked_trans[,2:8]

        all_trans_type <- rbind(true_swap,
                                linked_trans)

      } else {

        all_trans_type <- true_swap

      }

      all_trans_type["mig"] <- mig
      all_trans_type["param"] <- param
      all_trans_type["rep"] <- rep

    }

    all_trans_type["total_pairs_analysed"] <- nrow(W_80)
    all_trans_type["true_pairs_all"] <- nrow(trans_by_W80)

    all_trans_type <- add_observed_values(all_trans_type)

    all_phyloscanner_results <- rbind(all_phyloscanner_results, all_trans_type)

  }

  unlink("W_labels*")
  unlink("W_stats.csv")

}

#I will use the script below to create a function to sumarize the data



all_phyloscanner_results <- readRDS("all_phyloscanner_results_test2.RDS")
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


summarize_all_data(all_phyloscanner_results)


