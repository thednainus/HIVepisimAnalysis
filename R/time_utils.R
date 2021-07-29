# scripts to deal with time conversion



#' Title
#'
#' @param total_years total years simulated
#' @param last_sample_time date in which simulation would have taken place
#'    (let's say 2021)
#' @param sampleTimes sample times in days that will be converted to years
#' @param ini_sim date to used for the begining of simulations
#'
#' @return
#' @export
#'
#' @examples
days2years <- function(total_years, last_sample_time, sampleTimes, init_sim){

  time_sampled <- as.Date(sampleTimes, origin = init_sim)
  time_sampled <- decimal_date(time_sampled)

  return(time_sampled)
}




#' Title
#'
#' @param n_sample_times number of time to sample individuals
#'
#' @param art_init dataframe of ID and time on ART init
#' @param start_date start date for sampling times
#' @param end_date end date for sampling times
#' @param departures dataframe of ID and time of departures
#' @param diag_info dataframe of ID and time of diagnosis
#' @param tm transmission matrix
#' @param stages track history for stage of HIV infection of infected nodes
#' @param location location cab ne "region" or "global"
#'
#' @return
#' @export
#'
#' @examples
sampleIDs <- function(n_sample_times, start_date, end_date,
                      art_init, departures,
                      diag_info, stages, tm, location){

  #initially assume that we will have ids to sample
  ids2sample <- "yes"
  sampledIDs_list <- NULL

  #count for sample time
  count <- n_sample_times

  #get vector of IDs from location to sample from
  IDs <- get_id_by_location(tm = tm, location = location)

  #browser()
  while(ids2sample == "no" | count > 0){

      #sample a time within start_date and end_date
      sample_time <- runif(n = 1, min = start_date_dec, max = end_date_dec)

      # from IDs vector, get those IDs that have been diagnosed
      diag_ids <- diag_info[diag_info$IDs %in% IDs,]

      # then from diag_ids, get the ids that have been diagnosed before the sampled
      # time
      diagIDs_before_st <- diag_ids[diag_ids$time_decimal <= sample_time,]

      #from diagIDs_before_st, get ids that have departured before sampled time
      dep_before_st <- departured_ids(diagIDs_before_st, departures, sample_time)

      #from departured ids, get the ids that are active before sampled time
      active_ids <- diagIDs_before_st[!(diagIDs_before_st$IDs %in% dep_before_st$infID),]
      activeIDs_before_st <- active_ids[active_ids$time_decimal <= sample_time,]

      #from activeIDs_before_st, get ids that is on ART before sampled time
      ids_onART <- onART(activeIDs_before_st, art_init, sample_time)

      #from ids_onART, get the ones that are NOT on ART before sampled time
      ids_notOnArt <- activeIDs_before_st[!(activeIDs_before_st$IDs %in% ids_onART$IDs),]

      if(nrow(ids_notOnArt) == 0){
        ids2sample <- "no"
      }

      if(nrow(ids_notOnArt) > 0){

        #browser()

        results <- get_newids(ids_notOnArt$IDs)
        sid <- results[[1]]
        tmpids <- results[[2]]

        sampledIDs_list[[length(sampledIDs_list)+1]] <- set_sampledIDs(sid = sid,
                                                                       st = sample_time)
        print(paste("sampled id list", sampledIDs_list, sep = " "))

        index <- which(IDs == sid)
        IDs <- IDs[-index]

        count <- count - 1
        ids2sample <- "yes"

      }
  }
  return(sampledIDs_list)
}


#' Return the departures ID within the sampled time
#'
#' @param dig_ids_within_st dataframe of diagnosed ids and time
#' @param departures dataframe of departures
#' @param st sampled time
#'
#' @return
#' @export
#'
#' @examples
departured_ids <- function(dig_ids_within_st, departures, st){


  depIDs <- departures[departures$infID %in% dig_ids_within_st$IDs,]
  depIDs_within_st <- depIDs[depIDs$time_decimal <= st,]

  return(depIDs_within_st)

}


#' Title
#'
#' @param activeIDs_before_st
#' @param art_init
#' @param st sampled time
#'
#' @return
#' @export
#'
#' @examples
onART <- function(activeIDs_before_st, art_init, st){

  onARTdf <- art_init[art_init$IDs %in% activeIDs_before_st$IDs,]

  onART_IDs_within_st <- onARTdf[onARTdf$time_decimal <= st,]

  return(onART_IDs_within_st)

}


#' Return sampled id (sid) and vector of IDs to sample from
#'
#' @param ids ids to sample from
#'
#' @return
#' @export
#'
#' @examples
get_newids <- function(ids){


  if(length(ids) == 0){
    sid <- NULL
  }

  if(length(ids) == 1){
    sid <- ids
    tmpids <- NULL
  }

  if(length(ids) > 1) {
    sid <- sample(x = ids, size = 1)

    #sample another individual
    index <- which(ids == sid)
    tmpids <- ids[-index]
  }


  results <- list(sid, tmpids)
  return(results)
}







#' Get total number of individuals in network by region
#'
#' Get the total number of individuals in the final network by region of choice.
#'
#' @param tm transmission matrix
#' @param location location can be "region" or "global"
#'
#' @return
#' @export
#'
#' @examples
get_n_by_origin <- function(tm, location){

  if(location == "region"){
    new_tm <- subset(tm, infMigrant == 1 | infMigrant == 21)

  }

  if(location == "global"){
    new_tm <- subset(tm, infMigrant == 2 | infMigrant == 12)
  }

  el <- cbind(new_tm$inf, new_tm$sus)
  IDPOP <- unique(as.vector(el))

  n <- length(IDPOP)

  return(n)

}


#' Get the ID of individuals by location
#'
#' @param tm transmission matrix
#' @param location location can take the value of "region" or "global"
#'
#' @return
#' @export
#'
#' @examples
get_id_by_location <- function(tm, location){

  if(location == "region"){
    new_tm <- subset(tm, infMigrant == 1 | infMigrant == 21)

  }

  if(location == "global"){
    new_tm <- subset(tm, infMigrant == 2 | infMigrant == 12)
  }

  el <- cbind(new_tm$inf, new_tm$sus)
  IDPOP <- unique(as.vector(el))

  return(IDPOP)
}


#' Title
#'
#' @param sid The sampled ID
#' @param st The sampled time
#'
#' @return
#' @export
#'
#' @examples
set_sampledIDs <- function(sid, st){

  sampledID <- data.frame(sampled_ID = sid, sampled_time = st)

  return(sampledID)
}

