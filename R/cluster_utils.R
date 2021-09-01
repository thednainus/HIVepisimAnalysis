#functions to help generating pbs files to run analysis in the cluster

#' Function to split vector of with names of directories into chunks of
#' specific size
#'
#' @param dirs Vector with names of directories
#' @param size Size to divide vector.
#'
#' @details If we want to split a vector of 15 elements into chunks of 10 elements,
#'    the first object in the list will have 10 element, and the second object
#'    in the list will have 5 elements.
#'
#' @return
#' @export
split_dirs <- function(dirs, size){

  split_dirs <- split(dirs, ceiling(seq_along(dirs) / size))

  return(split_dirs)

}
