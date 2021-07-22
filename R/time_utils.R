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
  newids <- IDs[match(diagnosis$IDs, IDs)]
  newids <- newids[!is.na(newids)]

  sampledIDs_list <- NULL

  for(st in 1:length(sample_times)){

    #sample ID from diagnosed individuals only
    sid <- sample(x = newids, size = 1)

    # check whether sid was or ART at time of sampling
    if(any(art_init$IDs == sid) == TRUE){
      art_init_time <- art_init[art_init$IDs == sid,]

      if(nrow(art_init_time) == 1){
        #check whether individual was on ART at time of sampling
        if(sample_times[st] < art_init_time$time_decimal == TRUE){

          # then keep ID and sampled id
          sampledIDs_list[[st]] <- set_sampledIDs(sampled_ID = sid,
                                                  sampled_time = sample_times[st])

          #removed sampled id (sid) from vector of ids to sample from
          index <- which(newids == sid)
          newids <- newids[-index]
        }

        if(nrow(art_init_time) > 1){
          print("still have to do this part")
        }

        while(sample_times[st] < art_init_time$time_decimal == FALSE){

          #remove id from newids and sample another sid
          index <- which(newids == sid)
          tmpnewids <- newids
          tmpnewids <- tmpnewids[-index]

          sid <- sample(x = tmpnewids, size = 1)


        }



      }

    }




  }

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
