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
#' @param sample_times time in which individuals were sampled
#' @param art_init dataframe of ID and time on ART init
#' @param art_halt dataframe of ID and time on ART halt
#' @param art_reinit dataframe of ID and time on ART reinit
#' @param departures dataframe of ID and time of departures
#' @param diagnosis dataframe of ID and time of diagnosis
#' @param tm transmission matrix
#' @param stages track history for stage of HIV infection of infected nodes
#' @param location location cab ne "region" or "global"
#'
#' @return
#' @export
#'
#' @examples
sampleIDs <- function(sample_times, art_init, art_halt, art_reinit, departures,
                      diagnosis, stages, tm, location){

  IDs <- get_id_by_location(tm = tm, location = location)

  #from total IDPOP in TM
  #sample individuals that have been diagnosed
  #but when sampling, I still need to check whether individual was diagnosed
  #before or at the time of sampling
  newids <- IDs[match(diagnosis$IDs, IDs)]
  newids <- newids[!is.na(newids)]

  sampledIDs_list <- NULL

  for(st in 1:length(sample_times)){

    #sample ID from diagnosed individuals only
    sid <- sample(x = newids, size = 1)

    #this is just to create a loop
    #when keep_sid == TRUE
    keep_sid <- FALSE

    while(keep_sid == FALSE){

      if(any(art_init$IDs == sid) == TRUE){

        #check if individual was diagnosed before or at the time of sampling
        if(diagnosis$time_decimal[diagnosis$IDs == sid] <= sample_times[st]){

          #check if individual is on ART
          onART <- is_onART (sid = sid, st = sample_times[st], art_init = art_init)

          if(onART == "yes"){

            # check whether is on halt

            ART_or_halt <- is_onARThalt(sid = sid, st = sample_times[st], art_halt = art_halt)

            if(ART_or_halt == "onART"){

              #sampled another individual
              #TO DO
            }

            if(ART_or_halt == "haltART"){

              #check whether individual has initiated ART
              #sampled another individual
              #TO DO
            }

          }

          if(onART == "no"){

            # keep id if time of sampling was before individual started ART
            # then keep ID and sampled id
            sampledIDs_list[[st]] <- set_sampledIDs(sampled_ID = sid,
                                                    sampled_time = sample_times[st])

            #removed sampled id (sid) from vector of ids to sample from
            index <- which(newids == sid)
            newids <- newids[-index]

            keep_sid <- TRUE

          }
        }
      }
    }
  }
}


#' Check whether individual has re-initiated ART
#'
#' @param sid sampled id
#' @param st sampled time
#' @param art_reinit dataframe of art_reinit
#'
#' @return
#' @export
#'
#' @examples
is_ARTreinit <- function(sid, st, art_reinit){





}



#' Check whether individual has halt ART
#'
#' @param sid sampled id
#' @param st sampled time
#' @param art_halt dataframe with IDs and time that ART halt happened
#'
#' @return
#' @export
#'
#' @examples
is_onARThalt <- function(sid, st, art_halt){

  art_halt_time <- art_halt[art_halt$IDs == sid,]$time_decimal
  halt_index <- which(art_halt_time > st)

  if(length(halt_index) == 0){

    #in this condition individual is still on ART
    #this function will only be called if individual has initiated ART

    results <- "onART"


  }

  if(length(halt_index) > 0){

    halt_times <- art_halt_time[halt_index[1]]

    results <- "haltART"
  }

  return(results)

}

#' Check whether an individual is on ART
#'
#' @param sid sampled ID number
#' @param st sampled time
#' @param art_init dataframe for ART dates
#'
#' @return
#' @export
#'
#' @examples
is_onART <- function(sid, st, art_init){

  art_init_time <- art_init[art_init$IDs == sid,]

  # keep id if time of sampling was before individual started ART
  if(art_init_time$time_decimal < st){
    onART <- "yes"
  }else{
    onART <- "no"
  }

  return(onART)
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
#' @param sampled_ID The sampled ID
#' @param sampled_time The sampled time
#'
#' @return
#' @export
#'
#' @examples
set_sampledIDs <- function(sampled_ID, sampled_time){

  sampledID <- data.frame(sampled_ID = sid, sampled_time = sample_times[st])

  return(sampledID)
}
