#' Title Validate infector probability calculation
#'
#' This function saves several values to compare W calculated on the true trees
#'    and the true transmission matrix. For the return values see below.
#'
#' @param sim Simulation number to know the parameter values that network was
#'    simulated
#' @param run Run number
#' @param tm Transmission matrix per seed
#' @param W_true Infector probability calculated on the true trees
#' @param seed_ID Seed ID name. Seed is the node ID that started the
#'    transmissions; it is the initial infected node in the simulations.
#' @param MH maximum height in which infector probabilities were calculated
#' @param true_tree object of class phylo
#' @param prefix prefix to save results. If prefix = NULL, results will be saved
#'   on file W_stats.csv.
#'
#' @return a data frame of 1 row with the following values:
#'  1. Code = code indicating if W was calculated on true trees
#'  2. Maximum height = maximum height (MH) in which calculations were performed;
#'  3. Total tips = total number of tips in the phylogenetic tree;
#'  4. n_tips_region = total number of tips from region in the phylogenetic tree;
#'  5. n_tips_global = total number of tips from global in the phylogenetic tree;
#'  6. n_trans_W = total number of transmissions independent of infector probability
#'     values;
#'  7. n_trans_tm = total number of transmissions based on transmission matrix
#'     subset by MH;
#'  8. n_true_trans_W = total number of true transmission independent of infector
#'     probability (W) values. Here is the true transmission when compared to the
#'     transmission matrix;
#'  9. n_W80_all = number of transmission in which W >= 80%;
#'  10. n_W80_correctDonorRecipt = number of transmission in point 9 that is a
#'      true transmission (when comparing to the transmission matrix);
#'  11. n_W80_swapDonorRecipt = number of transmission in point 9 that is a
#'      true transmission (when comparing to the transmission matrix), but
#'      there was a swap between donor and recipient;
#'  12. expected_n_trans = expected number of transmission defined as the sum
#'      of all infector probability (W) values.
#'
#' @export
validadeW <- function(sim, run, tm, W_true, W_estimated, seed_ID, MH, true_tree, prefix = NULL){
  #get seed name
  IDnumber <- str_split(string = seed_ID, pattern = "_")[[1]][2]
  IDnumber <- paste("seed", IDnumber, sep = "_")

  #get tm by seed
  tm_seed <- tm[[IDnumber]]

  #add time in years
  tm_seed["time_years"] <- tm_seed["at"] * 1/365

  #subset tm to maximum height and to region only
  tm_mh <- subset(tm_seed, infOrigin == "region" & susOrigin == "region")

  # W on true trees ----
  W_stats <- W_manipulations(W_true, code = "W on true")


  # W on estimated trees ----
  # converting W_estimated on a dataframe
  W_on_estimated <- W_manipulations(W_estimated, code = "W on estimated")


  #converting tm_mh to the same collumn names as Wsub
  tm <- data.frame(donor_ID = tm_mh$inf, recip_ID = tm_mh$sus,
                   infectorProbability = 1, Code = "True transmission")


  # get transmission
  all_trans_Wtrue <- semi_join(W_stats$Wsub, tm, by = c("donor_ID", "recip_ID"))
  all_trans_West <- semi_join(W_on_estimated$Wsub, tm, by = c("donor_ID", "recip_ID"))




  if(!is.null(W_stats$W80)){
    # true transmissions: correct identification of donor and recipient
    W80true_correctDonorRecip <- semi_join(W_stats$W80, tm, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_stats$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    W80true_swapDonorRecip <- semi_join(W_stats$W80_trunc, tm, by = c("donor_ID", "recip_ID"))
  }

  if(!is.null(W_on_estimated$W80)){
    # true transmissions: correct identification of donor and recipient
    # based on estimated phylogenetic trees
    W80est_correctDonorRecip <- semi_join(W_on_estimated$W80, tm, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_on_estimated$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    # based on estimated phylogenetic trees
    W80est_swapDonorRecip <- semi_join(W_on_estimated$W80_trunc, tm, by = c("donor_ID", "recip_ID"))

  }



  #get names in phylogenetic trees
  region <- unlist(lapply(true_tree$tip.label, function(x) grepl("_1$", x) | grepl("_21", x)))
  global <- unlist(lapply(true_tree$tip.label, function(x) grepl("_2$", x) | grepl("_12", x)))


  all_data <- data.frame(sim,
                         run,
                         maximum_height = MH,
                         total_tips = length(true_tree$tip.label),
                         n_tips_region = sum(region),
                         n_tips_global = sum(global),
                         n_trans_Wtrue = nrow(W_stats$Wsub),
                         n_trans_West = nrow(W_on_estimated$Wsub),
                         n_trans_tm = nrow(tm),
                         n_true_trans_Wtrue = nrow(all_trans_Wtrue),
                         n_true_trans_West = nrow(all_trans_West),
                         n_W80true_all = ifelse(is.null(W_stats$W80), 0, nrow(W_stats$W80)),
                         n_W80true_correctDonorRecipt = ifelse(is.null(W_stats$W80),
                                                               0, nrow(W80true_correctDonorRecip)),
                         n_W80true_swapDonorRecipt = ifelse(is.null(W_stats$W80_trunc),
                                                            0, nrow(W80true_swapDonorRecip)),
                         n_W80est_all = ifelse(is.null(W_on_estimated$W80), 0, nrow(W_on_estimated$W80)),
                         n_W80est_correctDonorRecipt = ifelse(is.null(W_on_estimated$W80), 0,
                                                              nrow(W80est_correctDonorRecip)),
                         n_W80est_swapDonorRecipt = ifelse(is.null(W_on_estimated$W80_trunc), 0,
                                                           nrow(W80est_swapDonorRecip)),
                         expected_n_trans_Wtrue = sum(W_stats$W$infectorProbability),
                         expected_n_trans_West = sum(W_on_estimated$W$infectorProbability))

  if(is.null(prefix)){
    filename <- "W_stats.csv"
  } else {
    filename <- paste(prefix, "W_stats.csv", sep = "_")
  }

  write.table(all_data, file = filename, append = TRUE, sep = ",",
              row.names = FALSE, col.names = !file.exists(filename))

}


#' Summarize infector probability (W) for phylogenetic trees
#'
#' This function will summarize results for W calculated for a phylogenetic tree.
#'
#' @param sim Simulation number to know the parameter values that network was
#'    simulated.
#' @param tm Transmission matrix.
#' @param W1 Infector probability calculated on true trees (all region tips)
#' @param code string to identify if W was calculated on true or simulated trees.
#' @param prefix prefix to save results. If prefix = NULL, results will be saved
#'    as W.csv.
#' @param labels if labels = TRUE, labels of 0 or 1 will be created for ROC
#'    curves.
#'
#' @return a data frame of 1 row with the following values:
#'
#'   all_data <- data.frame(sim,
#'  1. Sim = simulation number (to get parameter values).
#'
#'
#'  2. Total tips = total number of tips in the phylogenetic tree.
#'
#'  3. n_tips_region = total number of tips from region in the phylogenetic tree.
#'
#'  4. n_tips_global = total number of tips from global in the phylogenetic tree.
#'
#'  5. n_trans_W_all = total number of transmission pairs in which W was calculated.
#'     Here we consider all infector probability values.Here W was calculated for
#'     the combination of tips in the phylogenetic tree.
#'
#'  6. n_trans_tm = total number of true transmissions based on the transmission
#'                  matrix (tm). Note that the tm will be subset to include only
#'                  IDs that are also observed in the phylogenetic tree.
#'
#'  7. n_true_trans_W = total number of true transmission independent of infector
#'     probability (W) values. Here is the true transmission when compared to the
#'     transmission matrix.
#'
#'  8. n_W80_all = number of transmission pairs in which W >= 80%.
#'
#'  9. n_W80_correctDonorRecipt = number of transmission pairs in point 8 that is a
#'      true transmission (when comparing to the transmission matrix).
#'
#'  10. n_W80_swapDonorRecipt = number of transmission pairs in point 8 that is a
#'      true transmission (when comparing to the transmission matrix), but
#'      there was a swap between donor and recipient.
#'
#'  11. expected_n_trans = expected number of transmission defined as the sum
#'                         of all infector probability (W) values.
#'
#' @export
summaryW <- function(sim, tm, W1, tree, code, prefix = NULL, labels = TRUE){

  #browser()

  #subset tm to region only
  tm_region <- subset(tm, infOrigin == "region" & susOrigin == "region")

  #keep_only_rows that tips are in phylogenetic tree
  rows_to_keep <- keep_row(df = tm_region, tree = tree)

  if(!is.null(rows_to_keep)){
    tm1 <- tm_region[rows_to_keep,]
  } else{
    #if there is no rows to keep, assign all values of tm1 to NA
    tm1 <- tm_region[1,]
    tm1[1,] <- NA
  }



  # W on true trees (all region tips: tree1) ----
  W_stats <- W_manipulations(W1, code = code)


  #converting tm to the same column names as Wsub
  # Wsub = data frame of W that was subset to contain only donor_ID and recip_ID
  # which we can compare to the transmission matrix
  tm_all1 <- data.frame(donor_ID = tm1$inf, recip_ID = tm1$sus,
                   infectorProbability = 1, Code = "True transmission")



  # get all transmissions that has a W calculated and also occurred as
  # true transmission in the transmission matrix
  all_trans_W1 <- semi_join(W_stats$Wsub, tm_all1, by = c("donor_ID", "recip_ID"))


  # W_stats$W80 is a subset of W that was higher than 80%
  if(!is.null(W_stats$W80)){
    # true transmissions: correct identification of donor and recipient
    W80true_correctDonorRecip <- semi_join(W_stats$W80, tm_all1, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_stats$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    W80true_swapDonorRecip <- semi_join(W_stats$W80_trunc, tm_all1, by = c("donor_ID", "recip_ID"))
  }



  #get names in phylogenetic trees
  region <- unlist(lapply(tree$tip.label, function(x) grepl("_1$", x) | grepl("_21", x)))
  global <- unlist(lapply(tree$tip.label, function(x) grepl("_2$", x) | grepl("_12", x)))



  all_data <- data.frame(sim,
                         total_tips = length(tree$tip.label),
                         n_tips_region_tree = sum(region),
                         n_tips_global_tree = sum(global),
                         n_trans_W_all = nrow(W_stats$Wsub),
                         n_trans_tm = nrow(tm1),
                         n_true_trans_W = nrow(all_trans_W1),
                         n_W80_all = ifelse(is.null(W_stats$W80), 0, nrow(W_stats$W80)),
                         n_W80_correctDonorRecipt = ifelse(is.null(W_stats$W80),
                                                               0, nrow(W80true_correctDonorRecip)),
                         n_W80_swapDonorRecipt = ifelse(is.null(W_stats$W80_trunc),
                                                            0, nrow(W80true_swapDonorRecip)),
                         expected_n_trans = sum(W_stats$W$infectorProbability))

  if(is.null(prefix)){
    filename <- "W_stats.csv"
  } else {
    filename <- paste(prefix, "W_stats.csv", sep = "_")
  }

  write.table(all_data, file = filename, append = FALSE, sep = ",",
              row.names = FALSE)

  if(labels == TRUE){

    ntips_all <- paste("all", length(tree$tip.label), sep = "")
    ntips_region <- paste("reg", sum(region), sep = "")
    ntips_global <- paste("glo", sum(global), sep = "")

    tip <- paste(ntips_all, ntips_region, ntips_global, sep = "_")

    if(is.null(prefix)){
      filename1 <- paste(paste("W_labels", tip, sep = "_"), "csv", sep = ".")
    } else {
      filename1 <- paste(paste(prefix, "W_labels", tip, sep = "_"), "csv", sep = ".")
    }

    # get transmission
    not_correct <- anti_join(W_stats$Wsub, tm_all1, by = c("donor_ID", "recip_ID"))
    not_correct["labels"] <- 0

    if(nrow(all_trans_W1) > 0){
      all_trans_W1["labels"] <- 1
    }


    all_labels <- rbind(all_trans_W1, not_correct)
    all_labels["total_tips"] <- length(tree$tip.label)
    all_labels["n_tips_region"] <- sum(region)
    all_labels["n_tips_global"] <- sum(global)

    #get only labels for cherries

    cherries <- get_tip_cherry(tree)
    cherries <- do.call(rbind, cherries)
    cherries <- as.data.frame(cherries)
    names(cherries) <- c("donor", "recip")

    cherries["donor_ID"] <- unlist(lapply(cherries$donor, function(x)
                                  str_split(string = x, pattern = "_")[[1]][1]))
    cherries["recip_ID"] <- unlist(lapply(cherries$recip, function(x)
                                  str_split(string = x, pattern = "_")[[1]][1]))

    cherries$donor_ID <- as.integer(cherries$donor_ID)
    cherries$recip_ID <- as.integer(cherries$recip_ID)


    cherries_truc <- data.frame(donor_ID = cherries$recip_ID, recip_ID = cherries$donor_ID)


    pairs12 <- semi_join(all_labels, cherries[3:4], by = c("donor_ID", "recip_ID"))
    pairs21 <- semi_join(all_labels, cherries_truc, by = c("donor_ID", "recip_ID"))

    pairs <- rbind(pairs12, pairs21)
    pairs$labels <- as.factor(pairs$labels)


    write.table(pairs, file = filename1, append = FALSE, sep = ",",
                row.names = FALSE)
  }

}



#' Summarize infector probability (W) for phylogenetic trees
#'
#' This function will summarize results for W calculated for a phylogenetic tree.
#' It will uses the information of sampling time. Pairs will be only considered
#' if transmission happened before sampling time.
#'
#' @param sim Simulation number to know the parameter values that network was
#'    simulated.
#' @param tm Transmission matrix.
#' @param W1 Infector probability calculated on true trees (all region tips)
#' @param ID ID name.
#' @param code String to identify if W was calculated on true or simulated trees.
#' @param st_df Dataframe for the IDs and its sampling times in decimal year.
#' @param init_sim_date Ititial date that simulation started in the format
#'    "1980-01-01".
#' @param prefix prefix to save results. If prefix = NULL, results will be saved
#'    as W.csv.
#' @param labels if labels = TRUE, labels of 0 or 1 will be created for ROC
#'    curves.
#'
#' @return a data frame of 1 row with the following values:
#'
#'   all_data <- data.frame(sim,
#'  1. Sim = simulation number (to get parameter values).
#'
#'
#'  2. Total tips = total number of tips in the phylogenetic tree.
#'
#'  3. n_tips_region = total number of tips from region in the phylogenetic tree.
#'
#'  4. n_tips_global = total number of tips from global in the phylogenetic tree.
#'
#'  5. n_trans_W_all = total number of transmission pairs in which W was calculated.
#'     Here we consider all infector probability values.Here W was calculated for
#'     the combination of tips in the phylogenetic tree.
#'
#'  6. n_trans_tm = total number of true transmissions based on the transmission
#'                  matrix (tm). Note that the tm will be subset to include only
#'                  IDs that are also observed in the phylogenetic tree.
#'
#'  7. n_true_trans_W = total number of true transmission independent of infector
#'     probability (W) values. Here is the true transmission when compared to the
#'     transmission matrix.
#'
#'  8. n_W80_all = number of transmission pairs in which W >= 80%.
#'
#'  9. n_W80_correctDonorRecipt = number of transmission pairs in point 8 that is a
#'      true transmission (when comparing to the transmission matrix).
#'
#'  10. n_W80_swapDonorRecipt = number of transmission pairs in point 8 that is a
#'      true transmission (when comparing to the transmission matrix), but
#'      there was a swap between donor and recipient.
#'
#'  11. expected_n_trans = expected number of transmission defined as the sum
#'                         of all infector probability (W) values.
#'
#' @export
summaryW2 <- function(sim, tm, W1, tree, code, st_df, init_sim_date,
                      prefix = NULL, labels = TRUE){

  #browser()

  #subset tm to region only
  tm_region <- subset(tm, infOrigin == "region" & susOrigin == "region")
  tm_region$year <- days2years(tm_region$at, init_date = init_sim_date)

  #keep_only_rows that tips are in phylogenetic tree
  rows_to_keep <- keep_row(df = tm_region, tree = tree)

  if(!is.null(rows_to_keep)){
    tm1 <- tm_region[rows_to_keep,]
  } else{
    #if there is no rows to keep, assign all values of tm1 to NA
    tm1 <- tm_region[1,]
    tm1[1,] <- NA
  }

  #pair of transmissions before sampling time
  tm_bst <- apply(st_df, 1, function(x) tm1[tm1$sus == x[1] & tm1$year <= x[2],])
  tm_bst <- do.call(rbind, tm_bst)

  #if(is.null(tm_bst)){
  #  browser()
  #  #if there is no rows to keep, assign all values of tm_bst to NA
  #  tm_bst <- tm_region[1,]
  #  tm_bst[1,] <- NA

  #}

  if(nrow(tm_bst) == 0){
    #if there is no rows to keep, assign all values of tm_bst to NA
    tm_bst <- tm_region[1,]
    tm_bst[1,] <- NA

  }



  # W on true trees (all region tips: tree1) ----
  W_stats <- W_manipulations(W1, code = code)


  #converting tm to the same column names as Wsub
  # Wsub = data frame of W that was subset to contain only donor_ID and recip_ID
  # which we can compare to the transmission matrix
  tm_all1 <- data.frame(donor_ID = tm_bst$inf, recip_ID = tm_bst$sus,
                        infectorProbability = 1, Code = "True transmission")



  # get all transmissions that has a W calculated and also occured as
  # true transmission in the transmission matrix
  all_trans_W1 <- semi_join(W_stats$Wsub, tm_all1, by = c("donor_ID", "recip_ID"))


  # W_stats$W80 is a subset of W that was higher than 80%
  if(!is.null(W_stats$W80)){
    # true transmissions: correct identification of donor and recipient
    W80true_correctDonorRecip <- semi_join(W_stats$W80, tm_all1, by = c("donor_ID", "recip_ID"))
  }
  if(!is.null(W_stats$W80_trunc)){
    # true transmission but incorrect identification of donor and recipient
    W80true_swapDonorRecip <- semi_join(W_stats$W80_trunc, tm_all1, by = c("donor_ID", "recip_ID"))
  }



  #get names in phylogenetic trees
  region <- unlist(lapply(tree$tip.label, function(x) grepl("_1$", x) | grepl("_21", x)))
  global <- unlist(lapply(tree$tip.label, function(x) grepl("_2$", x) | grepl("_12", x)))



  all_data <- data.frame(sim,
                         total_tips = length(tree$tip.label),
                         n_tips_region_tree = sum(region),
                         n_tips_global_tree = sum(global),
                         n_trans_W_all = nrow(W_stats$Wsub),
                         n_trans_tm = nrow(tm_bst),
                         n_true_trans_W = nrow(all_trans_W1),
                         n_W80_all = ifelse(is.null(W_stats$W80), 0, nrow(W_stats$W80)),
                         n_W80_correctDonorRecipt = ifelse(is.null(W_stats$W80),
                                                           0, nrow(W80true_correctDonorRecip)),
                         n_W80_swapDonorRecipt = ifelse(is.null(W_stats$W80_trunc),
                                                        0, nrow(W80true_swapDonorRecip)),
                         expected_n_trans = sum(W_stats$W$infectorProbability))

  if(is.null(prefix)){
    filename <- "W_stats_bst.csv"
  } else {
    filename <- paste(prefix, "W_stats_bst.csv", sep = "_")
  }

  write.table(all_data, file = filename, append = FALSE, sep = ",",
              row.names = FALSE)

  if(labels == TRUE){

    ntips_all <- paste("all", length(tree$tip.label), sep = "")
    ntips_region <- paste("reg", sum(region), sep = "")
    ntips_global <- paste("glo", sum(global), sep = "")

    tip <- paste(ntips_all, ntips_region, ntips_global, sep = "_")

    if(is.null(prefix)){
      filename1 <- paste(paste("W_labels_bst", tip, sep = "_"), "csv", sep = ".")
    } else {
      filename1 <- paste(paste(prefix, "W_labels_bst", tip, sep = "_"), "csv", sep = ".")
    }

    # get transmission
    not_correct <- anti_join(W_stats$Wsub, tm_all1, by = c("donor_ID", "recip_ID"))
    not_correct["labels"] <- 0

    if(nrow(all_trans_W1) > 0){
      all_trans_W1["labels"] <- 1
    }


    all_labels <- rbind(all_trans_W1, not_correct)
    all_labels["total_tips"] <- length(tree$tip.label)
    all_labels["n_tips_region"] <- sum(region)
    all_labels["n_tips_global"] <- sum(global)

    #get only labels for cherries

    cherries <- get_tip_cherry(tree)
    cherries <- do.call(rbind, cherries)
    cherries <- as.data.frame(cherries)
    names(cherries) <- c("donor", "recip")

    cherries["donor_ID"] <- unlist(lapply(cherries$donor, function(x)
      str_split(string = x, pattern = "_")[[1]][1]))
    cherries["recip_ID"] <- unlist(lapply(cherries$recip, function(x)
      str_split(string = x, pattern = "_")[[1]][1]))

    cherries$donor_ID <- as.integer(cherries$donor_ID)
    cherries$recip_ID <- as.integer(cherries$recip_ID)


    cherries_truc <- data.frame(donor_ID = cherries$recip_ID, recip_ID = cherries$donor_ID)


    pairs12 <- semi_join(all_labels, cherries[3:4], by = c("donor_ID", "recip_ID"))
    pairs21 <- semi_join(all_labels, cherries_truc, by = c("donor_ID", "recip_ID"))

    pairs <- rbind(pairs12, pairs21)
    pairs$labels <- as.factor(pairs$labels)


    write.table(pairs, file = filename1, append = FALSE, sep = ",",
                row.names = FALSE)
  }

}

W_manipulations <- function(W, code){

  #tm_mh$sus is the recipient
  #tm_mh$inf is the donor

  # W on true trees ----
  #convert W (infector probability matrix) to dataframe
  #browser()
  W1 <- as.data.frame(W)
  W1["donor_ID"] <- unlist(lapply(W$donor, function(x) str_split(x, pattern = "_")[[1]][1]))
  W1["recip_ID"] <- unlist(lapply(W$recip, function(x) str_split(x, pattern = "_")[[1]][1]))
  W1["Code"] <- code
  W1$donor_ID <- as.numeric(W1$donor_ID)
  W1$recip_ID <- as.numeric(W1$recip_ID)

  # subseting data frame W to contain only donor_ID and recip_ID
  # this is comparable to the transmission matrix
  Wsub <- data.frame(donor_ID = W1$donor_ID, recip_ID = W1$recip_ID,
                          infectorProbability = W1$infectorProbability, Code = W1$Code)
  Wsub$donor_ID <- as.integer(Wsub$donor_ID)
  Wsub$recip_ID <- as.integer(Wsub$recip_ID)

  # get only transmissions in which infector probability is >= 80%
  W80 <- Wsub[Wsub$infectorProbability >= 0.8,]
  # modify the order of recipient and donor just to check whether identification of
  # donor was not correct
  W80_trunc <- data.frame(donor_ID = W80$recip_ID, recip_ID = W80$donor_ID)


  return(list(W = W1, Wsub = Wsub, W80 = W80, W80_trunc = W80_trunc))


}


#' Get index of rows to keep in transmission matrix
#'
#' @param tm transmission matrix
#' @param tree object of class phylo
#'
#' @return Dataframe of rows to keep.
#' @export
keep_row <- function(df, tree){
  tip_IDs <- as.numeric(unlist(lapply(tree$tip.label, function(x) str_split(x, "_")[[1]][1])))
  rows_to_keep <- NULL

  #browser()
  while(length(tip_IDs) > 0){

    if(tip_IDs[1] %in% df$sus){
      index_sus <- which(tip_IDs[1] == df$sus)
      if(length(index_sus) == 1){
        if(df$inf[index_sus] %in% tip_IDs){
          rows_to_keep <- append(rows_to_keep, index_sus)
        }
      }
      if(length(index_sus) > 1 & any(df$inf[index_sus] %in% tip_IDs)){
        while(length(index_sus) > 0){
          if(df$inf[index_sus[1]] %in% tip_IDs){
            rows_to_keep <- append(rows_to_keep, index_sus[1])
            index_sus <- index_sus[-1]
          }else{index_sus <- index_sus[-1]}
        }
      }
    }

    if(tip_IDs[1] %in% df$inf){
      index_inf <- which(tip_IDs[1] == df$inf)
      if(length(index_inf) == 1){
        if(df$sus[index_inf] %in% tip_IDs){
          rows_to_keep <- append(rows_to_keep, index_inf)
        }
      }
      if(length(index_inf) > 1 & any(df$sus[index_inf] %in% tip_IDs)){
        while(length(index_inf) > 0){
          if(df$sus[index_inf[1]] %in% tip_IDs){
            rows_to_keep <- append(rows_to_keep, index_inf[1])
            index_inf <- index_inf[-1]
          }else{index_inf <- index_inf[-1]}
        }
      }
    }

    tip_IDs <- tip_IDs[-1]

  }

  rows_to_keep <- sort(unique(rows_to_keep))

  return(rows_to_keep)
}
