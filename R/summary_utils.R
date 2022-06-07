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
get_rates2 <- function(threshold, df_true){

  #browser()

  df_true["pred"] <- ifelse(df_true$infectorProbability >= threshold, "1", "0")
  df_true$pred <- as.factor(df_true$pred)

  df_true$labels <- as.character(df_true$labels)
  df_true$labels <- as.factor(df_true$labels)

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

    rates <- data.frame(threshold = threshold, TPR = TPR, FPR = FPR,
                        param = df_true$param, rep = df_true$rep,
                        perc = df_true$perc)

  }

  if(check_data_cm == "not_ok"){

    rates <- data.frame(threshold = threshold, TPR = NA, FPR = NA)

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


#' Calculate the log importance weight based on observed and simulated data
#'
#' @param diag_obs Frequency of observed incidence of HIV diagnosis in MSM
#'    from 1982 to 2020
#' @param diag_sim Frequency of estimated incidence of HIV diagnosis in MSM
#'    from 1982 to 2020
#'
#' @return
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
#' @return
#' @export
compute_log_importance_weight_incidence <- function ( incid_obs, incid_sim )
{
  #browser()
  incid_sim <- incid_sim$incid.pop1

  log_importance_weight <- sum( dpois( incid_sim, lambda = incid_obs , log = TRUE ) )

  return(log_importance_weight)
}
