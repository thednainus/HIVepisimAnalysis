#Summary utils: useful for plotting


#' Summarize data for statistics calculated for infector probability (W)
#'
#'
#' @param df dataframe containing data to summarize
#'
#' @return dataframe for median and quantiles for quantities of interest.
#'   see \code{\link{validadeW}}.
#' @export
summarize_data <- function(df){

  df <- df[c(1,5,9,18,19)]

  #convert df to long format
  df_long <- melt(df, id.vars = c("sim", "n_tips_region"))

  #calculate median by simulation
  median_by_sim <- aggregate(x = df_long$value ,                # Specify data column
                             by = list(df_long$sim, df_long$n_tips_region, df_long$variable),       # Specify group indicator
                             FUN = median)
  names(median_by_sim)[4] <- "median"

  #lower quantile
  lowerq_by_sim <- aggregate(x = df_long$value ,
                                by = list(df_long$sim, df_long$n_tips_region, df_long$variable),
                                FUN = quantile,
                                probs=0.025)
  names(lowerq_by_sim)[4] <- "lower"


  #upper quantile
  upperq_by_sim <- aggregate(x = df_long$value ,
                             by = list(df_long$sim, df_long$n_tips_region, df_long$variable),
                             FUN = quantile,
                             probs=0.975)
  names(upperq_by_sim)[4] <- "upper"

  all_data <- cbind(median_by_sim, lowerq_by_sim[4], upperq_by_sim[4])

  return(all_data)

}


#' Get epi data for simulation runs
#'
#' @param list_dirs list of directories containing data
#'
#' @return A list of dataframes for each type of simulations.
#'
#' @details This function will go over each simulation run saved as RDS format,
#'   convert it to a dataframe and assign as simulation number the number of the
#'   respective run. This function was designed when running more than 1
#'   simulation independently in the cluster, to merge all the runs together for
#'   plotting purposes.
#'
#' @export
get_epi_data <- function(list_dirs){

  #browser()

  for(n in 1:length(list_dirs)){

    if(n == 1){
      sim_all <- list()
    }

    type_name <- SplitPath(list_dirs[n])$filename
    list_runs <- dir(list_dirs[n], full.names = TRUE)

    count = 1
    for(r in 1:length(list_runs)){

      rds_sim <- list.files(list_runs[r], pattern = ".RDS", full.names = TRUE)
      #covert network sim into dataframe
      sim <- readRDS(rds_sim)
      #convert sim to data.frame
      sim_df <- as.data.frame(sim)
      sim_df$sim <- count
      sim_df$type <- type_name
      count <- count + 1

      if(r == 1){
        sim_dfs <- sim_df
      } else{sim_dfs <- rbind(sim_dfs, sim_df)}
    }
    sim_all[[n]] <- sim_dfs
  }
  return(sim_all)
}


#' Get rates to construct ROC curve
#'
#' Get false positive rates and true positive rate for ROC curves
#'
#' @param threshold list object with threshold values rangng from 0 to 1
#' @param df_true data.frame object for the infector probability and the true
#'    classification of transmission pairs
#'
#' @details The true data for df_true object will have 1 for a transmission pair
#'    that have occurred and 0 for a transmission pair that did not occur. This
#'    will be based by comparison with the transmission matrix.
#'
#' @return data.frame object for values of threshold, true positive rate (TPR)
#'    and false positive rate (FPR).
#' @export
get_rates <- function(threshold, df_true){

  #browser()

  df_true["pred"] <- ifelse(df_true$infectorProbability >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  #check whether the data can be used with the function caret:confusionMatrix

  check_data_cm <- check_data_cm(data = df_true$pred,
                                 reference = df_true$labels,
                                 positive = "1",
                                 mode = "sens_spec")

  if(check_data_cm == "ok"){

    #create confusion matrix
    cm <- confusionMatrix(df_true$pred, df_true$labels, positive = "1")
    #true positive rate
    TPR <- cm$byClass[[1]]
    #false positive rate
    FPR <- 1 - cm$byClass[[2]]

    rates <- data.frame(threshold = threshold, TPR = TPR, FPR = FPR)

  }

  if(check_data_cm == "not_ok"){

    rates <- data.frame(threshold = threshold, TPR = NA, FPR = NA)

  }

  return(rates)
}


#' Get rates to construct ROC curve
#'
#' Get false positive rates and true positive rate for ROC curves
#'
#' @param threshold list object with threshold values rangng from 0 to 1
#' @param df_true data.frame object for the infector probability and the true
#'    classification of transmission pairs
#'
#' @details The true data for df_true object will have 1 for a transmission pair
#'    that have occurred and 0 for a transmission pair that did not occur. This
#'    will be based by comparison with the transmission matrix.
#'
#' @return data.frame object for values of threshold, true positive rate (TPR)
#'    and false positive rate (FPR).
#' @export
get_rates2 <- function(threshold, df_true, byReplicate = FALSE){

  #browser()

  param <- unique(df_true$param)
  mig <- unique(df_true$mig)
  perc = unique(df_true$perc)
  code = unique(df_true$code)
  sampler = unique(df_true$sampler)
  mig = unique(df_true$mig)

  if(byReplicate == FALSE){
    rep = "merged"
  } else {
    rep = unique(df_true$rep)
  }


  df_true["pred"] <- ifelse(df_true$infectorProbability >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  df_true$labels <- as.character(df_true$labels)
  df_true$labels <- as.factor(df_true$labels)

  positives <- sum(df_true$labels == 1)
  negatives <- length(df_true$labels) - positives

  #check whether the data can be used with the function caret:confusionMatrix

  check_data_cm <- check_data_cm(data = df_true$pred,
                                 reference = df_true$labels,
                                 positive = "1",
                                 mode = "sens_spec")

  if(check_data_cm == "ok"){

    #create confusion matrix
    cm <- confusionMatrix(df_true$pred, df_true$labels, positive = "1")
    #true positive rate
    TPR <- cm$byClass[[1]]
    #false positive rate
    FPR <- 1 - cm$byClass[[2]]

    TN <- cm$table[1] #true negative
    FP <- cm$table[2] #false positive
    FN <- cm$table[3] #false negative
    TP <- cm$table[4] #true positive

    #browser()

    rates <- data.frame(threshold = threshold, TPR = TPR, FPR = FPR,
                        positives = positives,
                        negatives = negatives,
                        TN = TN,
                        FP = FP,
                        FN = FN,
                        TP = TP,
                        param = param, rep = rep,
                        perc = perc, code = code,
                        sampler = sampler, mig = mig)

  }

  if(check_data_cm == "not_ok"){


    rates <- data.frame(threshold = threshold, TPR = NA, FPR = NA,
                        param = param, rep = rep,
                        perc = perc, code = code,
                        sampler = sampler, mig = mig)

  }

  return(rates)
}


#' @export
check_data_cm <- function(data, reference, positive = "1", mode = "sens_spec"){

  # we assume here that the data is ok,
  #if it is not ok, the results object value will change in an if below
  results <- "ok"

  if (!(mode %in% c("sens_spec", "prec_recall", "everything")))
    stop("`mode` should be either 'sens_spec', 'prec_recall', or 'everything'")
  if (!is.factor(data) | !is.factor(reference)) {
    stop("`data` and `reference` should be factors with the same levels.",
         call. = FALSE)
  }
  if (!is.character(positive) & !is.null(positive))
    stop("positive argument must be character")
  if (length(levels(data)) > length(levels(reference))){
    #message("the data cannot have more levels than the reference")
    results <- "not_ok"
  }
  if (!any(levels(data) %in% levels(reference))) {
    #message("The data must contain some levels that overlap the reference.")
    results <- "not_ok"
  }
  if (!all(levels(data) %in% levels(reference))) {
    badLevel <- levels(data)[!levels(data) %in% levels(reference)]
    if (sum(table(data)[badLevel]) > 0) {
      #browser()
      #message("The data contain levels not found in the data.")
      results <- "not_ok"
    }
    #else {
    #  warning("The data contains levels not found in the data, but they are empty and will be dropped.")
    #  data <- factor(as.character(data))
    #}
  }
  if (any(levels(reference) != levels(data))) {
    #browser()
    #warning("HIVepisimAnalysis function:
    #        Levels are not in the same order for reference and data.
    #        Refactoring data to match.")
    data <- as.character(data)
    data <- factor(data, levels = levels(reference))
  }
  classLevels <- levels(data)
  numLevels <- length(classLevels)
  if (numLevels < 2){
    #message("there must be at least 2 factors levels in the data")
    results <- "not_ok"
  }


  return(results)
}


#' Get rates to construct ROC curve
#'
#' Get false positive rates and true positive rate for ROC curves
#'
#' @param threshold list object with threshold values rangng from 0 to 1
#' @param df_true data.frame object for the infector probability and the true
#'    classification of transmission pairs
#'
#' @details The true data for df_true object will have 1 for a transmission pair
#'    that have occurred and 0 for a transmission pair that did not occur. This
#'    will be based by comparison with the transmission matrix.
#'
#' @return data.frame object for values of threshold, true positive rate (TPR)
#'    and false positive rate (FPR).
#' @export
get_precision_recall <- function(threshold, df_true, byReplicate = FALSE){

  #browser()

  param <- unique(df_true$param)
  mig <- unique(df_true$mig)
  perc = unique(df_true$perc)
  code = unique(df_true$code)
  sampler = unique(df_true$sampler)
  mig = unique(df_true$mig)

  if(byReplicate == FALSE){
    rep = "merged"
  } else {
    rep = unique(df_true$rep)
  }


  df_true["pred"] <- ifelse(df_true$infectorProbability >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  df_true$labels <- as.character(df_true$labels)
  df_true$labels <- as.factor(df_true$labels)

  positives <- sum(df_true$labels == 1)
  negatives <- length(df_true$labels) - positives

  #check whether the data can be used with the function caret:confusionMatrix

  #check_data_cm <- check_data_cm(data = df_true$pred,
  #                               reference = df_true$labels,
  #                               positive = "1",
  #                               mode = "prec_recall")
  check_data_cm <- "ok"
  if(check_data_cm == "ok"){

    #create confusion matrix
    cm <- confusionMatrix(df_true$pred, df_true$labels,
                          positive = "1",
                          mode = "prec_recall")
    #true positive rate
    precision <- cm$byClass[[5]]
    #false positive rate
    recall <- cm$byClass[[6]]

    #browser()


    rates <- data.frame(threshold = threshold, precision = precision,
                        recall = recall,
                        positives = positives,
                        negatives = negatives,
                        param = param, rep = rep,
                        perc = perc, code = code,
                        sampler = sampler, mig = mig)

  }

  if(check_data_cm == "not_ok"){


    rates <- data.frame(threshold = threshold, precision = precision,
                        recall = recall,
                        positives = positives,
                        negatives = negatives,
                        param = param, rep = rep,
                        perc = perc, code = code,
                        sampler = sampler, mig = mig)

  }

  return(rates)
}

#' Calculate the log importance weight based on observed and simulated data
#'
#' @param diag_obs Frequency of observed incidence of HIV diagnosis in MSM
#'    from 1982 to 2020
#' @param diag_sim Frequency of estimated incidence of HIV diagnosis in MSM
#'    from 1982 to 2020
#'
#' @return Scalar for the log importance weight
#' @export
compute_log_importance_weight_newDx <- function ( diag_obs, diag_sim )
{
  #browser()
  diag_sim <- diag_sim$newDx_pop1
  #subset to remove the NA in observed data (data not available from 1981 and 2021)
  #diag_sim <- diag_sim[3:41]
  #diag_obs <- diag_obs[3:41]

  diag_sim <- diag_sim[6:41]
  diag_obs <- diag_obs[6:41]

  log_importance_weight <- sum( dpois( diag_sim, lambda = diag_obs , log = TRUE ) )

  return(log_importance_weight)
}

#' Calculate the log importance weight based on observed and simulated data
#' for incidence
#'
#' @param diag_obs Frequency of observed incidence in MSM from 1980 to 2020
#' @param diag_sim Frequency of estimated incidence in MSM from 1980 to 2020
#'
#' @return Scalar for the log importance weight
#' @export
compute_log_importance_weight_incidence <- function ( incid_obs, incid_sim )
{
  #browser()
  incid_sim <- incid_sim$incid.pop1

  log_importance_weight <- sum( dpois( incid_sim, lambda = incid_obs , log = TRUE ) )

  return(log_importance_weight)
}

#' Estimate difference between sampling time and time of infection
#'
#' Estimate the difference in months between time of infection of a susceptible
#' individual and time of sampling of that same individual. It assumes that a
#' month has 30 days.
#'
#' @param df1 Dataframe with information on time of sampling
#' @param df2 Transmission matrix to get information on when the individual
#'    got infected
#'
#' @return Dataframe with difference in months.
#' @export
get_difference <- function(df1, df2){
  #df1=dataframe of sampling times
  #df2=datafame of tm (transmission matrix)

  time_trans <- lapply(df1$sampled_ID, function(x) df2[df2$sus == x,][1:3])
  time_trans <- do.call(rbind, time_trans)

  time_trans["sampled_time"] <- df1$time_days
  time_trans["difference"] <- time_trans$sampled_time - time_trans$at
  time_trans["difference_months"] <- time_trans$difference/30

  return(time_trans)


}


#' Summarize results obtained with phyloscanner.
#'
#'
#' @param data
#'
#' @details Here we summarized the results obtained with phyloscanner.
#'    The summary is based on the summary results that phyloscanner returns
#'    in which it summarizes the "ancestry" type by window.
#'    We considered that two individuals are linked if proportion of values for
#'    ancestry is complex, multiTrans and trans are >= 0.50
#'    On the other hand, we considered that two individuals represent a true
#'    direction of transmission from host.1 to host.2 if the proportion of
#'    multiTrans and trans are >= 0.375. These cut-off values were based on
#'    Zhang et al.Clin Infect Dis. 2021 Jan 1; 72(1): 30–37.
#'    Note that to carry out this summary, we pre-removed all
#'    ancestry = noAncestry.
#'
#' @return Dataframe containing information of whether host.1 and host.2
#'    are linked and whether there direction of transmission is also correct.
#' @export
summarize_trans <- function(data){

  #browser()

  #data <- subset(data, ancestry != "noAncestry")

  #prop_linked <- sum(data$ancestry.tree.count)/(sum(data$both.exist)/nrow(data))
  prop_linked <- sum(data[data$ancestry != "noAncestry",]$ancestry.tree.count)/(sum(data$both.exist)/nrow(data))
  linked <- ifelse(prop_linked >= 0.5, "yes", "no")
  transmission <- data[data$ancestry == "trans" | data$ancestry == "multiTrans",]

  if(nrow(transmission) != 0){
    prop_trans <- sum(transmission$ancestry.tree.count)/(sum(transmission$both.exist)/nrow(transmission))
    #pairs were considered linked if prop_trans >= 0.375
    directness <- ifelse(prop_linked >= 0.5 & prop_trans >= 0.375, "yes", "no")


  } else {
    prop_trans <- 0
    directness <- "no"

  }

  prop_noAncestry <- sum(data[data$ancestry == "noAncestry",]$ancestry.tree.count)/(sum(data$both.exist)/nrow(data))

  pair_information <- tibble(prop_trans = prop_trans,
                             prop_linked = prop_linked,
                             prop_noAncestry = prop_noAncestry,
                             directness = directness,
                             linked = linked)

  return(pair_information)

}

#' Summarize results obtained with phyloscanner.
#'
#'
#' @param data
#'
#' @details Here we summarized the results obtained with phyloscanner.
#'    The summary is based on the summary results that phyloscanner returns
#'    in which it summarizes the "ancestry" type by window.
#'    We considered that two individuals are linked if proportion of values for
#'    ancestry is complex, multiTrans and trans are >= 0.50
#'    On the other hand, we considered that two individuals represent a true
#'    direction of transmission from host.1 to host.2 if the proportion of
#'    multiTrans, trans, and complex are >= 0.375.
#'    Here, we are interested to understand whether including complex transmissions
#'    would make it better to detect direction of transmission.
#'    Note that to carry out this summary, we pre-removed all
#'    ancestry = noAncestry.
#'
#' @return Dataframe containing information of whether host.1 and host.2
#'    are linked and whether there direction of transmission is also correct.
#' @export
summarize_trans_compl <- function(data){

  prop_linked <- sum(data$ancestry.tree.count)/(sum(data$both.exist)/nrow(data))
  linked <- ifelse(prop_linked >= 0.5, "yes", "no")
  transmission <- data[data$ancestry == "trans" | data$ancestry == "multiTrans",]

  if(nrow(transmission) != 0){
    prop_trans <- sum(transmission$ancestry.tree.count)/(sum(transmission$both.exist)/nrow(transmission))
    #pairs were considered linked if prop_trans >= 0.375
    direction <- ifelse(prop_trans >= 0.375, "yes", "no")


  } else {
    prop_trans <- 0
    direction <- "no"

  }

  pair_information <- tibble(prop_trans = prop_trans,
                             prop_linked = prop_linked,
                             direction = direction,
                             linked = linked)

  return(pair_information)

}



#' Check whether a pair of individuals represent a true transmission.
#'
#' This function will check whether a pair of individuals (host.1 and host.2)
#' represent a true transmission from host.1 to host.2 independent of phyloscanner
#' results.
#'
#' @param df1 Dataframe containing information of summarized results obtained
#'    with phyloscanner after running the function summarize_trans.
#' @param tm Dataframe containing information of the transmission matrix.
#'
#' @return Dataframe with a new column "trans" with only the individuals that
#'    represent a true transmission with correct direction of transmission.
#'    The value for column "trans" will be "true" (for true transmission).
#' @export
check_true_transmissions <- function(df1, tm){


  if("donor_ID" %in% names(df1)){

    true_transmissions <- tm[tm$inf == df1$donor_ID &
                               tm$sus == df1$recip_ID,]

  }else{
    true_transmissions <- tm[tm$inf == str_split(df1$host.1, "_")[[1]][2] &
                               tm$sus == str_split(df1$host.2, "_")[[1]][2],]
  }



  if(nrow(true_transmissions) == 0){
    df1["trans"] <- "false"

  }else if(nrow(true_transmissions) == 1){
    df1["trans"] <- "true"

  }else if(nrow(true_transmissions) > 1) {
    stop(print(nrow(true_transmissions)))
  }

  return(df1)

}


#' Check whether a pair of individuals represent a transmission pair.
#'
#' This function will check whether a pair of individuals (host.1 and host.2)
#' represent a transmission pair from host.1 to host.2 or host.2 to host.1
#' independent of phyloscanner results.
#'
#' @param df1 Dataframe containing information of summarized results obtained
#'    with phyloscanner after running the function summarize_trans.
#' @param tm Dataframe containing information of for the transmission matrix.
#'
#' @return Dataframe with a new column "real_pair" that take the value "true"
#'    (for a true transmission pair independent of directness,
#'    which means host.1 transmitted to host.2 or host.2 transmitted to host.1).
#'    It will also have another column "time_trans" for the time of transmission
#'    of a real pair. If a pair is not a real transmission pair, than the value
#'    for time_trans will be "NA".
#' @export
check_true_pair <- function(df1, tm){


  #this line was added to allow to sumarize data for consensus sequences
  if("donor_ID" %in% names(df1)){

    df1["inf"] <- df1$donor_ID
    df1["sus"] <- df1$recip_ID

  }else{
    df1["inf"] <- str_split(df1$host.1, "_")[[1]][2]
    df1["sus"] <- str_split(df1$host.2, "_")[[1]][2]
  }


  pairs <- tm[(tm$sus == df1$sus & tm$inf == df1$inf) |
                (tm$sus == df1$inf & tm$inf == df1$sus),]

  if(nrow(pairs) == 1){

    df1["real_pair"] <- "yes"
    df1["time_trans"] <- pairs$year

  }else if(nrow(pairs) == 0){

    df1["real_pair"] <- "no"
    df1["time_trans"] <- NA

  }else{
    #if more than one row appears in the dataframe pairs
    #something is wrong and should be further investigated.
    # This is because in my transmission model I only allowed infection from an
    #infected individual to a susceptible individual.
    stop(print(nrow))
  }

  if("donor_ID" %in% names(df1)){
    allpairs <- df1
  }else{
    allpairs <- df1[c(1:7,10:11)]
  }

  return(allpairs)

}


#' Check whether a pair of individuals represent a transmission pair.
#'
#' This function will check whether a pair of individuals (host.1 and host.2)
#' represent a transmission pair from host.1 to host.2 or host.2 to host.1
#' independent of phyloscanner results.
#'
#' @param df1 Dataframe containing information of summarized results obtained
#'    with phyloscanner after running the function summarize_trans.
#' @param tm Dataframe containing information of for the transmission matrix.
#'
#' @return Dataframe with a new column "real_pair" that take the value "true"
#'    (for a true transmission pair independent of directness,
#'    which means host.1 transmitted to host.2 or host.2 transmitted to host.1).
#'    It will also have another column "time_trans" for the time of transmission
#'    of a real pair. If a pair is not a real transmission pair, than the value
#'    for time_trans will be "NA".
#' @export
check_true_pair_directness <- function(df1, tm){


  df1["inf"] <- str_split(df1$host.1, "_")[[1]][2]
  df1["sus"] <- str_split(df1$host.2, "_")[[1]][2]

  pairs <- tm[(tm$sus == df1$sus & tm$inf == df1$inf),]

  if(nrow(pairs) == 1){

    df1["real_pair"] <- "yes"
    df1["time_trans"] <- pairs$year

  }else if(nrow(pairs) == 0){

    df1["real_pair"] <- "no"
    df1["time_trans"] <- NA

  }else{
    #if more than one row appears in the dataframe pairs
    #something is wrong and should be further investigated.
    # This is because in my transmission model I only allowed infection from an
    #infected individual to a susceptible individual.
    stop(print(nrow))
  }

  return(df1[c(1:7,10:11)])

}


#' Check whether a pair of individuals are linked by two or more individuals.
#'
#'
#' This function will check whether a pair of individuals (host.1 and host.2)
#' are linked because intermediary individuals were involved. For example, host.3
#' infected host.1 and host.2.
#'
#' @inheritParams check_true_transmissions
#' @param tm Dataframe of transmission matrix containing information of who
#'    transmitted to whom and at what time transmission happened.
#'
#' @return Dataframe with a new column "trans" with only the pair of individuals
#'    that are linked.
#'    The value for column "trans" will be "linked" or "not linked".
#' @export
check_linked_transmissions <- function(df1, tm){

  #browser()


  df1["inf"] <- str_split(df1$host.1, "_")[[1]][2]
  df1["sus"] <- str_split(df1$host.2, "_")[[1]][2]

  print(df1$sus)
  print(df1$inf)


  #tips_init <- unique(c(df1$host.1, df1$host.2))
  #tips_tm_init <- unlist(lapply(tips_init, function(x) str_split(x, pattern = "_")[[1]][2]))

  #tm_subset_list <- get_subset_tips(tips_tm_init, tm)

  transmissions <- get_chain(df1, tm)

  transmissions["inf"] <- df1$inf
  transmissions["sus"] <- df1$sus




  return(transmissions)

}

#' Get transmission chains
#'
#' @inheritParams  check_linked_transmissions
#' @param all_trans List of dataframes containing the transmission chains for
#'    each ID.
#' @param tm Transmission matrix containing the information of who transmitted
#'    whom.
#'
#' @details This function will go over each transmission to create a "transmission
#' chain from ID1 to ID5, for example. This will help indentify in which
#' circunstances phyloscanner indetified two IDs as a transmission pair,
#' when they were not one.
#' @return A dataframe with the transmission chain and the time of transmissions
#' @export
get_chain <- function(df1, tm){

  #convert dataframe to a graph using igraph
  df.g <- graph_from_data_frame(tm[,c(2:3,20)], directed = FALSE)

  #get shortest path between two nodes
  paths <- all_simple_paths(df.g, from = df1$inf, to = df1$sus)

  #get id names
  if(length(paths) == 1){
    chain <- as.vector(t(sapply(paths, as_ids)))

    #get years
    for(i in 1:(length(chain) - 1)){

      if(i == 1){
        year_trans <- tm[(tm$inf == chain[i] & tm$sus == chain[i+1]),]$year
        order <- paste(i, i+1, sep = "-")

        if(length(year_trans) == 0){
          year_trans <- tm[(tm$inf == chain[i+1] & tm$sus == chain[i]),]$year
          order <- paste(i+1, i, sep = "-")
        }
      } else{
        years <- tm[(tm$inf == chain[i] & tm$sus == chain[i+1]),]$year
        order2 <- paste(i, i+1, sep = "-")

        if(length(years) == 0){
          years <- tm[(tm$inf == chain[i+1] & tm$sus == chain[i]),]$year
          order2 <- paste(i+1, i, sep = "-")
        }

        year_trans <- c(year_trans, years)
        order <- c(order, order2)

      }

    }

  }

  if(length(paths) > 1){

    browser()

    print(paths)


  }

  if(length(paths) == 0){

    browser()


  }

  chain <- paste(chain, collapse = '-')
  year_trans <- paste(year_trans, collapse = '-')
  order <- paste(order, collapse = '_')


  chain_yearTrans <- data.frame(chain = chain, year_trans = year_trans,
                                order = order)

  return(chain_yearTrans)





  return(chain)


}



#' Get subset of susceptibles and infected individuals
#'
#' This function get the pairs of inf and sus individuals given a transmission
#' matrix.
#'
#' @inheritParams check_linked_transmissions
#' @param tips Vector of tip names (without the "ID")
#'
#' @return List of pairs of sus and inf individuals
#' @export
get_subset_tips <- function(tips, tm){

  tm_subset_list <- lapply(tips, function(x) tm[tm$sus == x | tm$inf == x ,])

  return(tm_subset_list)
}

#' Remove duplicated rows and order by date of transmission
#'
#' @param tm_subset_list List of dataframe of subset of transmission matrix
#'
#' @return Dataframe in which duplicated rows are removed and ordered by
#'    transmission date
#' @export
order_IDs_byDate <- function(tm_subset_list){


  #remove duplicated rows
  index <- NULL
  for(i in 1:length(tm_subset_list)){

    if(i < length(tm_subset_list) - 1){
      if(nrow(tm_subset_list[[i + 1]]) == 1){
        if(do.call(paste0, tm_subset_list[[i + 1]]) %in% do.call(paste0, tm_subset_list[[1]])){
          index <- c(index, i+1)
        }
        if(do.call(paste0, tm_subset_list[[i + 1]]) %in% do.call(paste0, tm_subset_list[[2]])){
          index <- c(index, i+1)
        }
      }
    }
  }
  index <- sort(unique(index))

  #remove elements from list

  if(!is.null(index)){
    tm_subset_final <- tm_subset_list[-index]
  }else{
    tm_subset_final <- tm_subset_list
  }

  # merge all into dataframe
  tm_subset_df <- do.call(rbind, tm_subset_final)
  #sort rows by year of transmission

  tm_subset_df <- tm_subset_df[order(tm_subset_df$sus, tm_subset_df$inf, tm_subset_df$year),]

  #remove duplicated rows
  tm_noDups <- tm_subset_df[!duplicated(tm_subset_df),]

  return(tm_noDups)
}


#' Summarize all results obtained with phyloscanner
#'
#' @param results_phylo Dataframe containing the values for number of transmission
#'    pairs analysed, number of true pairs within the total pairs analysed, number
#'    of true transmission pairs identified by phyloscanner, etc. by replicate,
#'    combination of parameter values and migration rate.
#'
#' @return Tibble object with the total number of pairs analysed, total number of
#'    true transmissions, total number of true positives as identified by phyloscanner,
#'    total number of false negatives by phyloscanner, total number of swaps
#'    (host.2 infected host.1, instead of host.1 infected host.2), and total
#'    number of false positives as identified by phyloscanner.
#'
#' @details This function is used in the R script summarize_phyloscanner_results_v2
#'    (check scripts in \href{https://github.com/thednainus/HIVepisimAnalysis/blob/main/Analyses/deep_sequencing/Result_analyses/summarize_phyloscanner_results_v2.R}{HIVepisimAnalysis GitHub})
#' @export
summarize_all_data <- function(results_phylo){

  total <- results_phylo %>%
    group_by(reps_groups) %>%
    summarize(total_pairs_analysed = total_pairs_analysed[[1]],
              true_pairs_all = sum(real_pair == "yes"))

  total_pairs_analysed <- sum(total$total_pairs_analysed)
  total_true_pairs_all <- sum(total$true_pairs_all)

  #total pairs that phyloscanner correctly identified the direction of transmissions
  true_phyloscanner <- subset(results_phylo, trans == "true" & direction == "yes")

  #total number of pairs that is a true transmission but phyloscanner identified
  #as not
  false_negative_phyloscanner <- subset(results_phylo, trans == "true" & direction == "no")

  #total number of pairs that phyloscanner identified as a transmission pair
  # but it is not
  false_trans_phyloscanner <- subset(results_phylo, (trans != "true" | is.na(trans)) & direction == "yes")


  #number of pairs that phyloscanner identified as correct direction of transmission
  #but it is a swap
  swap_phyloscanner <- subset(results_phylo, trans == "swap" & direction == "yes")
  if(nrow(swap_phyloscanner) == 0){
    swap_phyloscanner <- 0
  } else {
    swap_phyloscanner <- nrow(swap_phyloscanner)
  }

  summary_results <-data.frame(total_pairs_analysed = total_pairs_analysed,
                               total_true_pairs_all = total_true_pairs_all,
                               n_true_phyloscanner = nrow(true_phyloscanner),
                               n_false_negative_phyloscanner = nrow(false_negative_phyloscanner),
                               n_false_trans_phyloscanner = nrow(false_trans_phyloscanner),
                               swap_phyloscanner = swap_phyloscanner
  )

  return(summary_results)

}

#' Add observed values for phyloscanner results
#'
#' @param df1
#'
#' @return dataframe with the added observed values
#' @details   If observed value is TP,(true positives), it means that
#'    phyloscanner correctly identified the transmission pair. If observed value
#'    is FP (false positives), phyloscanner incorrectly identify the transmission
#'    pair. If observed value is FN (false negative) it means it is a true
#'    transmission pair, but phyloscanner was not able to identify it. If
#'    observed value is TN (true negative), it means it is not a transmission
#'    pair and phylosvanner correctly identified as not a transmission pair.
#' @export
add_observed_values <- function(df1){

  df1$observed <- ""

  df1$observed[df1$directness == "yes" &
                 df1$linked == "yes" &
                 df1$trans == "true"] <- "TP"

  df1$observed[df1$directness == "no" &
                 df1$linked == "yes" &
                 df1$trans == "true"] <- "FN"

  df1$observed[df1$directness == "no" &
                 df1$linked == "no" &
                 df1$trans == "true"] <- "FN"

  df1$observed[df1$directness == "no" &
                 df1$linked == "no" &
                 (is.na(df1$trans) | df1$trans == "swap")] <- "TN"

  df1$observed[df1$directness == "no" &
                 (df1$linked == "yes" | is.na(df1$linked)) &
                 (is.na(df1$trans) | df1$trans == "swap")] <- "TN"

  df1$observed[df1$directness == "yes" &
                 df1$linked == "yes" &
                 (is.na(df1$trans) | df1$trans == "swap")] <- "FP"

  return(df1)

}

#' Add observed values for phyloscanner results for direction
#'
#' @param df1
#'
#' @return dataframe with the added observed values
#' @details If observed value is TP (true positives), it means that
#'    phyloscanner correctly identified a transmission pair (independent of
#'    whether A infected B or B infected A). If observed value
#'    is FP (false positives), phyloscanner incorrectly identify a transmission
#'    pair. If observed value is FN (false negative), it means it is a true
#'    transmission pair, but phyloscanner was not able to identify it. If
#'    observed value is TN (true negative), it means it is not a transmission
#'    pair and phyloscanner correctly identified it as not a transmission pair.
#' @export
add_observed_values_direction <- function(df1){

  df1$observed <- ""

  df1$observed[df1$real_pair == "yes" &
                 df1$linked == "yes"] <- "TP"

  df1$observed[df1$real_pair == "no" &
                 df1$linked == "yes"] <- "FP"

  df1$observed[df1$real_pair == "yes" &
                 df1$linked == "no"] <- "FN"

  df1$observed[df1$real_pair == "no" &
                 df1$linked == "no"] <- "TN"



  return(df1)

}


#' Add observed values for phyloscanner results for directness
#'
#' @param df1
#'
#' @return dataframe with the added observed values
#' @details If observed value is TP (true positives), it means that
#'    phyloscanner correctly identified a transmission pair with correct direction
#'    of transmission(A infected B or B infected A). If observed value
#'    is FP (false positives), phyloscanner incorrectly identify a transmission
#'    pair. If observed value is FN (false negative), it means it is a true
#'    transmission pair (with correct direction of transmission from A to B or
#'    B to A), but phyloscanner was not able to identify it. If
#'    observed value is TN (true negative), it means it is not a transmission
#'    pair and phyloscanner correctly identified it as not a transmission pair
#'    (with correct direction of transmission from A to B or B to A).
#' @export
add_observed_values_directness <- function(df1){

  df1$observed <- ""

  df1$observed[df1$real_pair == "yes" &
                 df1$directness == "yes"] <- "TP"

  df1$observed[df1$real_pair == "no" &
                 df1$directness == "yes"] <- "FP"

  df1$observed[df1$real_pair == "yes" &
                 df1$directness == "no"] <- "FN"

  df1$observed[df1$real_pair == "no" &
                 df1$directness == "no"] <- "TN"



  return(df1)

}

#' Function to add the sampled times for transmission pairs
#'
#' @param df1 dataframe containing information of transmission pairs.
#' @param sampled_times dataframe containing information of sampled times
#'    for IDs in df1.
#'
#' @return dataframe with added information of the sampled times for donor
#'    (infected individual) and recipient (susceptible individual).
#' @details This function will adataframe composed of a single row and check
#'    whether it represents a transmission pair (independent of A transmitted
#'    to B or B transmitted to A). If it is a transmision pair will add to the
#'    dataframe the sampled times of infected and susceptible individuals.
#'    If not a transmission pair, than the value of NA will be added to the
#'    dataframe.
#' @export
add_sampled_times <- function(df1, sampled_times){

  if(df1$trans == "true"){

    #host.1 is the donor
    df1["st_donor_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.1, "_")[[1]][2],]$sampled_time
    df1["st_recip_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.2, "_")[[1]][2],]$sampled_time


  } else if(df1$trans == "false" & df1$real_pair == "yes"){

    #host.2 in this case is the donor
    df1["st_donor_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.2, "_")[[1]][2],]$sampled_time
    df1["st_recip_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.1, "_")[[1]][2],]$sampled_time


  } else if(df1$real_pair == "no"){
    #even if it is not a real pair, there is a sampled time associated to the IDs.
    #But because it is not a real pair, we won't have a time of transmission.
    df1["st_donor_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.1, "_")[[1]][2],]$sampled_time
    df1["st_recip_ID"] <- sampled_times[sampled_times$sampled_ID == str_split(df1$host.2, "_")[[1]][2],]$sampled_time

  }


  return(df1)



}


#' Who infected whom
#'
#' @param df1 dataframe containing information of the pairs that a TP (true positives),
#'    and the correct direction of transmission.
#'
#' @details This function will return the absolute numbers for the total number
#'    of pairs that are TP and from those pairs the total number in which
#'    phyloscanner also correctly estimated who infected whom.
#'
#' @return dataframe containg the absolute number of the total pairs that are
#'    TP and from those the total number of pairs in which phyloscanner correctly
#'    estimated who infected whom.
#' @export
wiw <- function(df1){

  #browser()

  TP <- subset(df1, observed == "TP")

  correct_direction <- subset(TP, real_pair == "yes" & directness == "yes")

  total_correct <- nrow(correct_direction)
  total_pairs <- nrow(TP)

  absolute_numbers <- data.frame(total_TP = total_pairs,
                                 total_correct = total_correct)

  return(absolute_numbers)


}

