# scripts to deal with sampling

#' Convert days to years
#'
#' This function converts days to years taking into consideration the start date
#' of simulations.
#'
#' @param sampleTimes Vector for sample times in days
#' @param init_date The date in the format of year-month-day to convert days to
#'    years
#'
#' @return Vector of sample times in decimal years
#' @export
#'
#' @examples
#' sample_times <- c(8760, 10526, 14600)
#' init_date <- "1981-01-01"
#' days2years(sample_times, init_date)
days2years <- function(sampleTimes, init_date){

  time_sampled <- base::as.Date(sampleTimes, origin = init_date)
  time_sampled <- lubridate::decimal_date(time_sampled)

  return(time_sampled)
}

#' Sample IDs of individuals diagnosed, active and not on ART
#'
#' Function to sample IDs based on sampled time and metadata. This function will
#' sample individuals that are diagnosed, active and not on ART at the time of
#' sampling.
#'
#' @param perc Percentage of the population to be sampled: Numbers should be between
#'    0 and 1.
#' @param start_date Start date for sampling times in decimal year.
#' @param end_date End date for sampling times in decimal year.
#' @param art_init Dataframe of IDs and time each ID initiated antiretroviral
#'    treatment (ART). Time must be in decimal year.
#' @param departure Dataframe of IDs and their time of departure. Time must be
#'    in decimal year.
#' @param diag_info Dataframe of IDs and their time of diagnosis. Time must be
#'    in decimal year.
#' @param tm Transmission matrix of who transmitted to whom.
#' @param location Character of location at time of transmission to get IDs
#'    to sample from. Location can take the value of "region" or "global".
#'
#' @details This function will get all IDs in transmission matrix (tm) from
#'    location. From those IDs, it selected all IDs that have been diagnosed
#'    during a simulation run. From all diagnosed IDs, it will select all IDs
#'    that have been diagnosed before sampled time. Then, from this new list of
#'    IDs, it select all IDs that have been active before time of sampling.
#'    Finally, from active and diagnosed IDs before sampled times, it will get
#'    all IDs that never initiated ART before sample times. From those IDs not
#'    on ART, the function will random sample a ID to be associated to the
#'    sampled time. This function is to mimic the idea of sampling individuals
#'    on HIV studies.
#'
#' @return Dataframe of IDs and their sampled times
#' @export
sampleIDs <- function(perc, start_date, end_date,
                      art_init, departure,
                      diag_info, tm, location){

  #browser()

  #initially assume that we will have ids to sample
  ids2sample <- "yes"
  sampledIDs_list <- NULL

  #vector of sids already sampled
  sampled_sids <- NULL


  #get vector of IDs from location to sample from
  IDs <- get_id_by_location(tm = tm, location = location)

  # from IDs vector, get those IDs that have been diagnosed
  diag_ids <- diag_info[diag_info$IDs %in% IDs,]

  #estimate the total number of individuals to sample based on percentage
  #of total diagnosed and active individuals
  #count for sample time
  count <- est_sampleSize(perc, start_date, end_date,
                          art_init, departure,
                          diag_info, tm, location, onART = TRUE)


  #browser()
  while(ids2sample == "no" | count > 0){

      #sample a time within start_date and end_date
      sample_time <- runif(n = 1, min = start_date, max = end_date)

      # then from diag_ids, get the ids that have been diagnosed before the sampled
      # time
      diagIDs_before_st <- diag_ids[diag_ids$time_decimal <= sample_time,]

      #from diagIDs_before_st, get ids that have departed before sampled time
      dep_before_st <- selectIDs(diagIDs_before_st$IDs, departure, sample_time)

      #from departed ids, get the ids that are active before sampled time
      active_ids_before_st <- diagIDs_before_st[!(diagIDs_before_st$IDs %in% dep_before_st$infID),]
      #activeIDs_before_st <- active_ids[active_ids$time_decimal <= sample_time,]

      #from activeIDs_before_st, get ids that is on ART before sampled time
      ids_onART <- selectIDs(activeIDs_before_st$IDs, art_init, sample_time)


      #from ids_onART, get the ones that are NOT on ART before sampled time
      ids_notOnArt <- activeIDs_before_st[!(activeIDs_before_st$IDs %in% ids_onART$IDs),]

      if(nrow(ids_notOnArt) == 0){
        ids2sample <- "no"
      }

      if(nrow(ids_notOnArt) > 0){

        # this will make to sure to remove the sids that hev already been sampled
        if(is.null(sampled_sids)){
          sid <- get_newids(ids_notOnArt$IDs)
          sampled_sids <- c(sampled_sids, sid)
        }else{
          ids_notOnArt <- ids_notOnArt[! ids_notOnArt$IDs %in% sampled_sids,]
          sid <- get_newids(ids_notOnArt$IDs)
          sampled_sids <- c(sampled_sids, sid)
        }

        #tmpids <- results[[2]]

        sampledIDs_list[[length(sampledIDs_list)+1]] <- set_sampledIDs(sid = sid,
                                                                       st = sample_time)
        #print(paste("sampled id list", sampledIDs_list, sep = " "))

        #index <- which(IDs == sid)
        #IDs <- IDs[-index]

        count <- count - 1
        ids2sample <- "yes"

      }
  }
  sampledIDsdf <- do.call(rbind, sampledIDs_list)
  return(sampledIDsdf)
}


#' Sample IDs of individuals who have been diagnosed, are active and on or not on ART
#'
#' Function to sample IDs based on sampled time and metadata. This function will
#' sample individuals that are diagnosed and active (individuals can or cannot be on ART)
#' at the time of sampling.
#'
#' @param perc Percentage of the population to be sampled: Numbers should be between
#'    0 and 1.
#' @param start_date Start date for sampling times in decimal year.
#' @param end_date End date for sampling times in decimal year.
#' @param art_init Dataframe of IDs and time each ID initiated antiretroviral
#'    treatment (ART). Time must be in decimal year.
#' @param departure Dataframe of IDs and their time of departure. Time must be
#'    in decimal year.
#' @param diag_info Dataframe of IDs and their time of diagnosis. Time must be
#'    in decimal year.
#' @param tm Transmission matrix of who transmitted to whom.
#' @param location Character of location at time of transmission to get IDs
#'    to sample from. Location can take the value of "region" or "global".
#'
#' @details This function will get all IDs in transmission matrix (tm) from
#'    location. From those IDs, it selected all IDs that have been diagnosed
#'    during a simulation run. From all diagnosed IDs, it will select all IDs
#'    that have been diagnosed before sampled time. Then, from this new list of
#'    IDs, it select all IDs that have been active before time of sampling.
#'    Finally, from active and diagnosed IDs before sampled times, it will get
#'    all IDs that never initiated ART before sample times. From those IDs not
#'    on ART, the function will random sample a ID to be associated to the
#'    sampled time. This function is to mimic the idea of sampling individuals
#'    on HIV studies.
#'
#' @return Dataframe of IDs and their sampled times
#' @export
sampleIDs2 <- function(perc, start_date, end_date,
                      art_init, departure,
                      diag_info, tm, location){

  #browser()

  #initially assume that we will have ids to sample
  ids2sample <- "yes"
  sampledIDs_list <- NULL

  #vector of sids already sampled
  sampled_sids <- NULL

  #get vector of IDs from location to sample from
  IDs <- get_id_by_location(tm = tm, location = location)

  # from IDs vector, get those IDs that have been diagnosed
  diag_ids <- diag_info[diag_info$IDs %in% IDs,]

  #get percentage of total number of individuals for each population
  #estimate the total number of individuals to sample based on percentage
  #of total diagnosed individuals => start_date
  diag_ids_start_date <- diag_ids[diag_ids$time_decimal >= start_date,]
  n_sample_times <- round(perc * nrow(diag_ids_start_date), 0)

  #browser()
  while(ids2sample == "no" | count > 0){

    #sample a time within start_date and end_date
    sample_time <- runif(n = 1, min = start_date, max = end_date)

    # then from diag_ids, get the ids that have been diagnosed before the sampled
    # time
    diagIDs_before_st <- diag_ids[diag_ids$time_decimal <= sample_time,]

    #from diagIDs_before_st, get ids that have departed before sampled time
    dep_before_st <- selectIDs(diagIDs_before_st$IDs, departure, sample_time)

    #from departed ids, get the ids that are active before sampled time
    active_ids_before_st <- diagIDs_before_st[!(diagIDs_before_st$IDs %in% dep_before_st$infID),]
    #activeIDs_before_st <- active_ids[active_ids$time_decimal <= sample_time,]

    #from activeIDs_before_st, get ids that is on ART before sampled time
    #ids_onART <- selectIDs(activeIDs_before_st$IDs, art_init, sample_time)


    #from ids_onART, get the ones that are NOT on ART before sampled time
    #ids_notOnArt <- activeIDs_before_st[!(activeIDs_before_st$IDs %in% ids_onART$IDs),]

    if(nrow(activeIDs_before_st) == 0){
      ids2sample <- "no"
    }

    if(nrow(activeIDs_before_st) > 0){

      #browser()

      # this will make to sure to remove the sids that hev already been sampled
      if(is.null(sampled_sids)){
        sid <- get_newids(activeIDs_before_st$IDs)
        sampled_sids <- c(sampled_sids, sid)
      }else{
        activeIDs_before_st <- activeIDs_before_st[! activeIDs_before_st$IDs %in% sampled_sids,]
        sid <- get_newids(activeIDs_before_st$IDs)
        sampled_sids <- c(sampled_sids, sid)
      }


      #tmpids <- results[[2]]

      sampledIDs_list[[length(sampledIDs_list)+1]] <- set_sampledIDs(sid = sid,
                                                                     st = sample_time)
      #print(paste("sampled id list", sampledIDs_list, sep = " "))

      #index <- which(IDs == sid)
      #IDs <- IDs[-index]

      count <- count - 1
      ids2sample <- "yes"

    }
  }
  sampledIDsdf <- do.call(rbind, sampledIDs_list)
  return(sampledIDsdf)
}


#' Get IDs before sampled time.
#'
#' @param IDs Vector of IDs
#' @param metadata Dataframe of specific metadata to get IDs before sampled time.
#' @param st Sampled time in decimal year.
#'
#' @details This function will get a vector of IDs and get the IDs that
#'    fit the condition of existing before the sampled time.
#'
#' @return Dataframe of IDs and the time of interest in decimal year.
#' @export
selectIDs <- function(IDs, metadata, st){

  selectedIDs <- metadata[metadata$IDs %in% IDs,]

  selectedIDs_before_st <- selectedIDs[selectedIDs$time_decimal <= st,]

  return(selectedIDs_before_st)

}


#' Get a sampled ID (sid) and vector of IDs to sample from
#'
#' @param ids Vector of ids to sample from
#'
#' @return One element for sampled ID (ID).
#'
#' @details This function will get a vector of IDs to sample from and randomly
#'    sample an ID. If there is only one element in ids, then sid will be equal
#'    to ids.
#'
#' @seealso \code{\link[base]{sample}} for more information on random samples.
#'
#' @export
#'
get_newids <- function(ids){

  if(length(ids) == 1){
    sid <- ids
  }

  if(length(ids) > 1) {
    sid <- sample(x = ids, size = 1)

  }

  return(sid)
}


#' Get total number of individuals in network by region
#'
#' Get the total number of individuals in the final network by region of choice.
#'
#' @inheritParams sampleIDs
#'
#' @return Integer of the total number of individuals in location of choice.
#'    Note that the function will return the total number of individual at location
#'    of choice at time of seroconversion.
#'
#' @details Note that location will be those at the time of seroconvertion.
#' @export
get_n_by_origin <- function(tm, location){

  if(location == "region"){
    new_tm <- subset(tm, susMigrant == 1 | susMigrant == 21)

  }

  if(location == "global"){
    new_tm <- subset(tm, susMigrant == 2 | susMigrant == 12)
  }

  #el <- cbind(new_tm$inf, new_tm$sus)
  IDPOP <- new_tm$sus

  n <- length(IDPOP)

  return(n)

}


#' Get ID of individuals by location
#'
#' This function will get the ID of infected individuals from specific location
#' (region or global) at the time of seroconversion.
#'
#' @inheritParams sampleIDs
#'
#' @return Vector of IDs from location of choice.
#'
#' @details Note that location will be those at the time of seroconvertion.
#' @export
get_id_by_location <- function(tm, location){

  if(location == "region"){
    new_tm <- subset(tm, susMigrant == 1 | susMigrant == 21)

  }

  if(location == "global"){
    new_tm <- subset(tm, susMigrant == 2 | susMigrant == 12)
  }

  #el <- cbind(new_tm$inf, new_tm$sus)
  IDPOP <- new_tm$sus

  return(IDPOP)
}


#' Set a dataframe of sampled IDs (sid) and their sampled times.
#'
#' @param sid The sampled ID
#' @inheritParams selectIDs
#'
#' @return Return dataframe of sampled IDs and their sampled times.
#' @export
set_sampledIDs <- function(sid, st){

  sampledID <- data.frame(sampled_ID = sid, sampled_time = st)

  return(sampledID)
}

#' Estimation sample size based on percentage
#'
#' @inheritParams sampleIDs
#'
#' @details This function will estimate the sample size based on a percentage
#' (0 to 1) and the end of the simulation (end_date). To carry out our analyses,
#' we used a sampled to sample time of sequencing just to mimic real HIV studies.
#' However, we can only sample diagnosed and active (and on ART) inidviduals,
#' it is hard to estimate the exact number of individuals based on the sampling
#' times that will be randomly sampled.
#' @return
#' @export
est_sampleSize <- function(perc, start_date, end_date, art_init,
                           departure, diag_info, tm, location = "region",
                           onART = FALSE){


  #get vector of IDs from location to sample from
  IDs <- get_id_by_location(tm = tm, location = location)

  # from IDs vector, get those IDs that have been diagnosed
  diag_ids <- diag_info[diag_info$IDs %in% IDs,]

  #get percentage of total number of individuals for each population
  #estimate the total number of individuals to sample based on percentage
  #of total diagnosed individuals => start_date
  diag_ids_bf_end_date <- diag_ids[diag_ids$time_decimal <= end_date,]

  #from diagIDs_before_st, get ids that have departed before sampled time
  dep_before_st <- selectIDs(diagIDs_before_st$IDs, departure, end_date)

  #from departed ids, get the ids that are active before sampled time
  active_ids <- diagIDs_before_st[!(diagIDs_before_st$IDs %in% dep_before_st$infID),]

  if(onART == TRUE){
    #from activeIDs_before_st, get ids that is on ART before sampled time
    ids_onART <- selectIDs(active_ids$IDs, art_init, end_date)
    #from ids_onART, get the ones that are NOT on ART before sampled time
    ids_notOnArt <- active_ids[!(active_ids$IDs %in% ids_onART$IDs),]

    n_sample_times <- round(perc * nrow(ids_notOnArt), 0)

  }else{
    n_sample_times <- round(perc * nrow(active_ids), 0)
  }

  return(n_sample_times)

}




