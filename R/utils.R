#' @title Get transmat infection tree for each seed or for specified seeds
#'
#' @description Separate a \code{transmat} object into different seedling
#'       event.
#'
#' @param x An object of class \code{"transmat"}, the output from
#'        \code{\link{get_transmat}}.
#' @param vertex.exit.times  optional numeric vector providing the time of
#'        departure of vertices, to be used to scale the lengths of branches
#'        reaching to the tips. Index position on vector corresponds to network
#'        id. NA indicates no departure, so branch will extend to the end of the
#'        tree.
#' @param by_areas by_areas = "all" if interested in all seeds, or by_areas = "region"
#'        if interested in seeds for region only.
#' @param max_value max value for comparison to select the seed values. It will
#'        select only seed values less than max_value for seeds from region only.
#'
#' @details
#' Converts a \code{\link{transmat}} object containing information about the
#' history of a simulated infection into a separated transmission matrix by
#' seedling event.
#'
#' This function was created using as base the function
#' \code{\link{as.phylo.transmat}}
#'
#'
#' @export
#'

get.transmat.phylo <- function(x, vertex.exit.times, by_areas = "all", max_value = NULL) {

 by_areas <-  by_areas
 max_value <- max_value
  # if not named properly, assume inf, sus at
  if (!all(c("inf", "sus", "at") %in% names(x))) {
    warning("input does not have appropriate column names for transmat,
            assuming first 3 should be 'inf','sus','at'")
    names(x) <- c("inf", "sus", "at")
  }
  tm <- x
  if (missing(vertex.exit.times)) {
    vertex.exit.times <- NULL
  }
  # find roots (infectors that never appear as sus)
  if(by_areas == "all"){
    v <- setdiff(unique(tm$inf), unique(tm$sus))
  }
  if(by_areas == "region"){
    v <- setdiff(unique(tm$inf), unique(tm$sus))
    v <- v[v < max_value]
  }
  if (length(v) > 1) {
    message("found multiple trees, returning a list of ", length(v),
            "phylo objects")
    # need to extract the portions of the edgelist and call seperately
    subtm <- lapply(v, function(v_sub) {
      # walk down the list to find elements below v_sub
      sub_rows <- which(tm$inf == v_sub)
      toFind <- v_sub
      while (length(toFind) > 0) {
        i <- toFind[1]
        sub_rows <- unique(c(sub_rows, which(tm$inf == i)))
        toFind <- c(toFind[-1], tm$sus[which(tm$inf == i)])
      }
      # call as.phylo on the subset of the edgelist
      get.transmat.phylo(tm[sub_rows, , drop = FALSE],
                         vertex.exit.times = vertex.exit.times,
                         by_areas = by_areas,
                         max_value = max_value)
      if (exists("subtm")){
        subtm <- rbind(subtm, tm[sub_rows, , drop = FALSE])
      }else{
        subtm <- tm[sub_rows, , drop = FALSE]
      }

    })
    names(subtm) <- paste("seed", v, sep = "_")
    return(subtm)
  }
}





#' Create transmission matrix csv
#'
#' @description Create a transmission matrix file to be used with the
#'    VirusTreeSimulator. VirusTreeSimulator requires that seeds are included
#'    in the file.
#'
#' @inheritParams create_sample_csv
#' @param time_tr Time of transmission for seeds. The other time of transmission
#'    should be provided in the tm dataframe.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as inf.csv
#'
#' @return Write data to csv file.
#' @export
create_inf_csv <- function(tm, time_tr, prefix = NULL){

  seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

  if(length(time_tr) == 1){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)),
                            rep(1980, length(seed_names)))
  }else if(length(time_tr)!= 1 & length(time_tr) == length(seed_names)){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)), time_tr)
  }else{
    stop("`time_tr` should be of length 1 or length of number of seeds")
  }


  colnames(seed_idtr) <- c("IDREC", "IDTR", "TIME_TR")

  #inf_sus <- data.frame(tm$sus, tm$inf, (tm$at * (1/365)))
  inf_sus <- data.frame(tm$sus, tm$inf, tm$year)
  colnames(inf_sus) <- c("IDREC", "IDTR", "TIME_TR")

  all_data <- rbind(seed_idtr, inf_sus)

  if(is.null(prefix)){
    filename <- "inf.csv"
  } else {
    filename <- paste(prefix, "inf.csv", sep = "_")
  }

  write.csv(x = all_data, file = filename, row.names = FALSE)
}


#' Create transmission matrix csv
#'
#' @description Create a transmission matrix file to be used with the
#'    VirusTreeSimulator. VirusTreeSimulator requires that seeds are included
#'    in the file.
#'
#' @inheritParams create_sample_csv
#' @param time_tr Time of transmission for seeds. The other time of transmission
#'    should be provided in the tm dataframe.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as inf.csv
#'
#' @return Write data to csv file.
#' @export
create_inf_csv2 <- function(tm, time_tr, prefix = NULL){

  seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

  if(length(time_tr) == 1){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)),
                            rep(0, length(seed_names)))
  }else if(length(time_tr)!= 1 & length(time_tr) == length(seed_names)){
    seed_idtr <- data.frame(seed_names, rep(NA, length(seed_names)), time_tr)
  }else{
    stop("`time_tr` should be of length 1 or length of number of seeds")
  }


  colnames(seed_idtr) <- c("IDREC", "IDTR", "TIME_TR")

  inf_sus <- data.frame(tm$sus, tm$inf, (tm$at * (1/365)))
  #inf_sus <- data.frame(tm$sus, tm$inf, tm$year)
  colnames(inf_sus) <- c("IDREC", "IDTR", "TIME_TR")

  all_data <- rbind(seed_idtr, inf_sus)

  if(is.null(prefix)){
    filename <- "inf.csv"
  } else {
    filename <- paste(prefix, "inf.csv", sep = "_")
  }

  write.csv(x = all_data, file = filename, row.names = FALSE)
}

#' Get sample csv file
#'
#' In this function the sample file for VirusTreeSimulator program
#' will contain all individuals in the transmision matrix
#'
#' @param tm Transmission matrix as returned using the function
#'    \link[EpiModel]{get_transmat}
#' @param time_seq Vector for sample time. If one value is provided,
#'    it will be replicated to the size of samples in the transmission matrix (tm).
#' @param seq_count The number of sequences per sample in the transmission matrix.
#'    If one value is provided, it will be replicated to the size of samples in
#'    the transmission matrix (tm).
#' @param prefix Text for prefix to use when saving filename.
#'
#' @details
#' If a prefix is not provided, csv file will be saved as sample.csv
#'
#' @return Write data to csv file.
#' @export
#'
#' @examples
#' To Do
create_sample_csv <- function(tm, time_seq, seq_count, prefix = NULL){

  IDPOP <- union(tm$inf, tm$sus)

  if(length(time_seq) == 1){
    TIME_SEQ <- rep(time_seq, length(IDPOP))
  } else if (length(time_seq) == length(IDPOP)){
    TIME_SEQ <- time_seq
  }else if (length(time_seq) != length(IDPOP)){
    stop("`time_seq` should be of length 1 or length of samples in the transmission matrix")
  }


  if(length(seq_count) == 1){
    SEQ_COUNT <- rep(seq_count, length(IDPOP))
  } else if (length(seq_count) == length(IDPOP)){
    SEQ_COUNT <- seq_count
  }else if(length(seq_count) != length(IDPOP)){
    stop(" `seq_count` should be of length 1 or length of samples in the transmission matrix")
  }


  all_data <- data.frame(IDPOP, TIME_SEQ, SEQ_COUNT)

  if(is.null(prefix)){
    filename <- "sample.csv"
  } else {
    filename <- paste(prefix, "sample.csv", sep = "_")
  }


  write.csv(x = all_data, file = filename, row.names = FALSE)

}


#' Get sample csv file
#'
#' In this function, the sample file for VirusTreeSimulator program
#' will contain only the sampled individuals (to mimic a real HIV study)
#'
#' @param ids Vector of samples IDs.
#' @inheritParams create_sample_csv
#'
#' @details
#' If a prefix is not provided, csv file will be saved as sample.csv
#'
#' @return Write data to csv file.
#' @export
create_sample_csv2 <- function(ids, time_seq, seq_count, prefix = NULL){

  IDPOP <- ids

  if(length(time_seq) == 1){
    TIME_SEQ <- rep(time_seq, length(IDPOP))
  } else if (length(time_seq) == length(IDPOP)){
    TIME_SEQ <- time_seq
  }else if (length(time_seq) != length(IDPOP)){
    stop("`time_seq` should be of length 1 or length of samples in
         the transmission matrix")
  }


  if(length(seq_count) == 1){
    SEQ_COUNT <- rep(seq_count, length(IDPOP))
  } else if (length(seq_count) == length(IDPOP)){
    SEQ_COUNT <- seq_count
  }else if(length(seq_count) != length(IDPOP)){
    stop(" `seq_count` should be of length 1 or length of samples
         in the transmission matrix")
  }


  all_data <- data.frame(IDPOP, TIME_SEQ, SEQ_COUNT)

  if(is.null(prefix)){
    filename <- "sample.csv"
  } else {
    filename <- paste(prefix, "sample.csv", sep = "_")
  }


  write.csv(x = all_data, file = filename, row.names = FALSE)

}


#' Get tsv for times
#'
#' In this function, the times file for CoaTrans program
#' will contain only the sampled individuals (to mimic a real HIV study)
#'
#' @param ids Vector of samples IDs.
#' @inheritParams create_sample_csv
#'
#' @details
#' If a prefix is not provided, tab delimited file file will be saved as times.tsv
#'
#' @return Write data to tsv file.
#' @export
create_times_tsv <- function(ids, time_seq, seq_count, prefix = NULL){

  IDPOP <- ids

  if(length(time_seq) == 1){
    TIME_SEQ <- rep(time_seq, length(IDPOP))
  } else if (length(time_seq) == length(IDPOP)){
    TIME_SEQ <- time_seq
  }else if (length(time_seq) != length(IDPOP)){
    stop("`time_seq` should be of length 1 or length of samples in
         the transmission matrix")
  }

  all_data <- data.frame(IDPOP, TIME_SEQ)

  if(seq_count > 1){
    rows <- c(1:nrow(all_data))
    times <- seq_count
    all_data <- all_data[rep(rows, times),]
  }

  if(is.null(prefix)){
    filename <- "times.tsv"
  } else {
    filename <- paste(prefix, "times.tsv", sep = "_")
  }


  write.table(x = all_data, file = filename, sep = "\t",
              row.names = FALSE, col.names = FALSE)

}


#' Create transmission matrix tsv
#'
#' @description Create a transmission matrix file to be used with the
#'    CoaTrans. CoaTrans. requires that seeds are included
#'    in the file.
#'
#' @inheritParams create_sample_csv
#' @param time_tr Time of transmission for seeds. The other time of transmission
#'    should be provided in the tm dataframe.
#'
#' @details
#' If a prefix is not provided, tsv file will be saved as transmissions.tsv
#'
#' @return Write data to tsv file.
#' @export
create_transmissions_tsv <- function(tm, time_tr, prefix = NULL){

  seed_names <- setdiff(unique(tm$inf), unique(tm$sus))

  if(length(time_tr) == 1){
    seed_idtr <- data.frame(rep("None", length(seed_names)), seed_names,
                            rep(1980, length(seed_names)))
  }else if(length(time_tr)!= 1 & length(time_tr) == length(seed_names)){
    seed_idtr <- data.frame(rep("None", length(seed_names)), seed_names, time_tr)
  }else{
    stop("`time_tr` should be of length 1 or length of number of seeds")
  }


  colnames(seed_idtr) <- c("IDTR", "IDREC", "TIME_TR")

  #inf_sus <- data.frame(tm$sus, tm$inf, (tm$at * (1/365)))
  sus_inf <- data.frame(tm$inf, tm$sus, tm$year)
  colnames(sus_inf) <- c("IDTR", "IDREC", "TIME_TR")

  all_data <- rbind(seed_idtr, sus_inf)

  if(is.null(prefix)){
    filename <- "transmissions.tsv"
  } else {
    filename <- paste(prefix, "transmissions.tsv", sep = "_")
  }

  write.table(x = all_data, file = filename, quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
}



#' @export
cd4s <- function(stage){

  if (stage==0) return(1e3)
  if (stage==1) return(750)
  if (stage==2) return(400)
  if (stage==3) return(300)
  if (stage==4) return(100)
}


#' Get CD4s by node IDs
#'
#' @param df_actives dataframe of IDs that were active at the end of a simulation
#'    and stages of HIV infection.
#' @param df_departures dataframe of IDs that departed the network before end of
#'    simulation and stages of HIV infection
#'
#' @return Vector of CD4 count.
#' @export
#'
get_cd4s <- function(IDPOP, df_actives, df_departures){

  df_actives["cd4s"] <- unlist(lapply(df_actives$stage_inf, function(x) cd4s(x)))
  #cd4s of active nodes at the end of simulation
  active_cd4s <- df_actives$cd4s
  active_cd4s_IDs <- setNames(active_cd4s, df_actives$infID)

  df_departures["cd4s"] <- unlist(lapply(df_departures$stage_dep, function(x) cd4s(x)))
  dep_cd4s <- df_departures$cd4s
  dep_cd4s_IDs <- setNames(dep_cd4s, df_departures$infID)

  index1 <- match(IDPOP, names(active_cd4s_IDs))
  cd4s1 <- active_cd4s_IDs[index1]
  index2 <- match(IDPOP, names(dep_cd4s_IDs))
  cd4s2 <- dep_cd4s_IDs[index2]

  index_tmp <- index1
  index_tmp[is.na(index_tmp)] <- index2[!is.na(index2)]
  #get cds4
  cd4s_tmp <- cd4s1

  cd4s_tmp[is.na(cd4s_tmp)] <- cd4s2[!is.na(cd4s2)]
  names(cd4s_tmp)[is.na(names(cd4s_tmp))] <- names(cd4s2)[!is.na(names(cd4s2))]

  return(cd4s_tmp)

}


#' Get CD4s by node IDs at time of sampling
#'
#' @param sampled_info dataframe of sampled IDs and sampled times
#' @param stages dataframe that tracks the stage of HIV infection for each ID.
#'
#' @details This function will get the stage of HIV infection at time of sampling.
#'
#' @return Vector of CD4 count.
#' @export
#'
get_cd4s_sampling <- function(sampled_info, stages){


  stages_per_ID <- apply(sampled_info, MARGIN = 1,
                         function(x) tail(subset(stages, infIDs == as.numeric(x["sampled_ID"]) &
                                              time_decimal <= as.numeric(x["sampled_time"])), n = 1))
  stages_per_ID <- do.call(rbind, stages_per_ID)


  stages_per_ID["cd4"] <- unlist(lapply(stages_per_ID$HIVstages, function(x) cd4s(x)))

  cd4_count <- stages_per_ID$cd4
  cd4_count <- setNames(cd4_count, stages_per_ID$infIDs)


  return(cd4_count)

}


#' Get which IDs are at early stage of HIV infection in the past 6 months
#'
#' This function mimics the idea of the recency test. It selects the IDs of
#' individuals with recent HIV infection (in the past 6 months)
#'
#' @inheritParams get_cd4s_sampling
#'
#' @return Vector of IDs.
#' @export
recency_test <- function(sampled_info, stages){

  stages_per_ID <- apply(sampled_info, MARGIN = 1,
                         function(x) tail(subset(stages, infIDs == as.numeric(x["sampled_ID"]) &
                                                   time_decimal <= as.numeric(x["sampled_time"])), n = 1))
  #this is the HIV stage of infection of the ID at time of sampling
  stages_per_ID <- do.call(rbind, stages_per_ID)
  #order dataframe to get the difference and time considered on early HIV stage
  #of infection
  #select only individuals on early HIV stage
  stages0 <- stages_per_ID[stages_per_ID$HIVstages == 0,]
  stages0_ordered <- stages0[order(stages0$infIDs),]

  #get the sampled times and IDs of individuals in stage0
  sampled_info_stage0 <- sampled_info[sampled_info$sampled_ID %in% stages0_ordered$infIDs,]

  #order sampled IDs and sampled time
  sampled_info_stage0_ordered <- sampled_info_stage0[order(sampled_info_stage0$sampled_ID),]

  sampled_info_stage0_ordered["difference"] <- sampled_info_stage0_ordered$sampled_time - stages0_ordered$time_decimal

  #return IDs in which the difference is equal or less than 0.5 (representing
  # half of the year or approx. 6 months)
  recent_HIV_infection <- sampled_info_stage0_ordered[sampled_info_stage0_ordered$difference <= 0.5,]

  return(recent_HIV_infection$tip_name)

}



