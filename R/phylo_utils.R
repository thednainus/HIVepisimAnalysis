#' Get tip names of phylogenetic tree based on transmission matrix
#'
#' @inheritParams get.transmat.phylo
#' @param format If format = "origin" return tip in the form of ID_global or ID_region.
#'        If format = "migrant" return tip in the form of ID_1, ID_2, ID_21, ID_12
#'
#' @return List of tip names. See details for codes on the tip names related
#'    to origin or global.
#'
#' @details
#' This code was build using the algorithm as in \code{\link{as.phylo.transmat}}.
#' The difference is that it will return only the tip names of phylogenetic tree.
#' If using format = "origin", tip names will be in the form of ID_global or
#' ID_region. If using format = "migrant", tip names will be in the form of
#' ID_1, ID_21, ID_2, ID_12.
#' Migrant values can take the value of 1 if vertex is from region,
#' 2 if vertex is from global,
#' 21 if vertex is from global that migrated to region and
#' 12 if vertex is from region that migrated to global.
#' It will return the origin at time of seroconversion.
#'
#'
#' This script will specifically work with the package HIVepisim.
#'
#' @export
get_tipNames <- function(x, format = "migrant", by_areas = "region",
                          max_value = NULL) {
  format <- format
  by_areas <-  by_areas
  max_value <- max_value

  # if not named properly, assume inf, sus at
  if (!all(c("inf", "sus", "at") %in% names(x))) {
    warning("input does not have appropriate column names for transmat,
            assuming first 3 should be 'inf','sus','at'")
    names(x) <- c("inf", "sus", "at")
  }
  tm <- x

  # find roots (infectors that never appear as sus)
  if(by_areas == "all"){
    v <- setdiff(unique(tm$inf), unique(tm$sus))
  }
  if(by_areas == "region"){
    v <- setdiff(unique(tm$inf), unique(tm$sus))
    v <- v[v < max_value]
  }

  if(length(v) >= 1){
    el <- cbind(tm$inf, tm$sus)
    origNodes <- unique(as.vector(el))

    if (format == "origin"){
      el_ori <- cbind(at = tm$at,
                      inf = tm$inf,
                      sus = tm$sus,
                      inf_ori = paste(tm$inf, tm$infOrigin, sep = "_"),
                      sus_ori = paste(tm$sus, tm$susOrigin, sep = "_"))
    }

    if (format == "migrant"){
      el_ori <- cbind(at = tm$at,
                      inf = tm$inf,
                      sus = tm$sus,
                      inf_ori = paste(tm$inf, tm$infMigrant, sep = "_"),
                      sus_ori = paste(tm$sus, tm$susMigrant, sep = "_"))
    }


    # get tip names in the form of inf_infOri or sus_susOri
    #browser()
    tipNamesOri <- unique(as.vector(el_ori[,c(4:5)]))
    # get the ID name of tips without the origin. This will allow to determine duplication
    # when a tip was region and changed to global (or vice-versa, for example)
    # or when using the different values for the "migrant" option (that can be 1
    # for region, 2 for global, 21 indicating a global that migrated to region and
    # 12 indicating a region that migrated to global)
    tip_names <- unlist(lapply(tipNamesOri, function(x) str_split(string = x, pattern = "_")[[1]][1]))
    # get the duplicated
    # add unique function to get only one value per duplicated
    # fromLast = TRUE gets the first value of the duplicated
    dup_tip <- unique(tip_names[duplicated(tip_names, fromLast = TRUE)])

    if(length(dup_tip) > 0){
      tip_origin_names <- get_newTip_names(dup_tip, el_ori)
      final_tip_names <- replace_dups(dup_tip, origNodes, tipNamesOri, tip_origin_names, tip_names)
    } else {
      final_tip_names <- tipNamesOri
    }

    return(final_tip_names)
  }
}



#' Get the time of duplicates at seroconversion.
#'
#' @description This is a support function for the main function
#'    \code{\link{toPhylo_transmatOrigin}}.
#'
#' @param dup_rows matrix of list of matrices of transmat
#'    that have duplicated vertex IDs.
#' @param get_time vector of oldest time
#'
#' @details During the simulations, vertex can migrate from region to global
#' or global to region. When adding the information of migrant or origin to the
#' tip name, we will observe duplicated on vertex IDs because it can be in region
#' at one time point but later in global, for example.
#' This function will get the duplicated rows in the transmission matrix
#' in which a duplicated vertex is observed. It will then get the time observed
#' for each duplicated vertex based on get_time. Basically it will return the
#' origin at time of seroconversion.
#'
#' @export
time_atSeroconversion <- function(dup_rows, get_time){
  #browser()
  dup_rows <- as.data.frame(dup_rows)
  row_int <- dup_rows[dup_rows$at == get_time,]

  if(nrow(row_int) > 1){
    browser()
  }

  return(row_int)
}

#' Get tip names based on duplicated rows.
#'
#' @description This is a support function for the main function
#'    \code{\link{toPhylo_transmatOrigin}}.
#'
#' @param dup_tip Vector of duplicated tips
#' @param el_ori Matrix of edge list and information on origin or migration
#'
#' @return A vector of the duplicated vertex in the form of ID_origin or ID_migrant.
#'    Based on the time of duplicated rows at seroconversion
#'     (see \code{\link{time_atSeroconversion}}),
#'    it will return the tip name that should be in the phylogenetic tree.
#'    I used the origin or migrant information observed at the time of
#'    seroconversion.
#'
#' @export
get_newTip_names <- function(dup_tip, el_ori){
  #browser()
  if(length(dup_tip) == 1){
    # get duplicated rows in the el_ori matrix
    dup_rows <- el_ori[el_ori[,2] == dup_tip | el_ori[,3] == dup_tip,]
    get_time <- head(sort(as.numeric(dup_rows[,1])), 1)
    row_int <- time_atSeroconversion(dup_rows, get_time)

    tip_origin_names <- name_of_tips(dup_tip, row_int)

  } else if (length(dup_tip) > 1 ){
    # get duplicated rows in the el_ori matrix
    dup_rows <- lapply(dup_tip, function(x) el_ori[el_ori[,2] == x | el_ori[,3] == x,])
    # if duplicated tip appears in more the one row, get the oldest time
    get_time <- lapply(dup_rows, function(x) head(sort(as.numeric(x[,1])), 1))
    #get the row of interest
    row_int <- mapply(time_atSeroconversion, dup_rows, get_time, SIMPLIFY = FALSE)

    # here get unique values because we can still observe at a same time point a
    # transmission by a same node (individual)
    tip_origin_names <- unique(unlist(unname(mapply(name_of_tips, dup_tip, row_int))))
  }

  return(tip_origin_names)

}

#' @title Get tip names.
#'
#' @description This is a support function for the main function
#'    \code{\link{toPhylo_transmatOrigin}}.
#'
#' Function that gets the name of the tip in the form of ID_origin or ID_migrant
#' because we can observe duplicates (when in the form ID_origin or ID_migrant),
#' we need to remove those from the phylogenetic tree.
#'
#' @inheritParams get_newTip_names
#' @param row_int data.frame of the row of interest. The row in the transmission
#' matrix that have the oldest time step based on duplicated value.
#' See \code{\link{time_atSeroconversion}} and \code{\link{get_newTip_names}}.
#'
#' @return vector of ID_origin or ID_migrant.
#'
#' @export
name_of_tips <- function(dup_tip, row_int){


  if(dup_tip %in% row_int$sus){
    infs_values <- row_int$sus_ori
  }else if(dup_tip %in% row_int$inf){
    infs_values <- row_int$inf_ori
  }

  return(infs_values)

}

#' @title Replace duplicates.
#'
#' @description This is a support function for the main function
#'    \code{\link{toPhylo_transmatOrigin}}.
#'
#' Replace duplicated in the form o ID_origin or ID_migrant by the correct tip
#' name.
#'
#' @inheritParams get_newTip_names
#'
#' @param origNodes Vector of vertex IDs.
#' @param tipNamesOri Vector in the form of ID_origin or ID_migrant. This vector
#'    contain the duplicated vertex when considering only the ID.
#' @param tip_names_toInsert Vector of the correct value of tip name to use in
#'    the phylogenetic tree in the form of ID_origin or ID_migrant
#' @param tip_names Vector of tip names in the form of ID only. This vector will
#'    contain the duplicated ID or IDs.
#'
#' @return Vector of tip names in the form of ID_origin or ID_migrant to use
#'    as tip.labels in the phylogenetic tree.
#'
#' @export
replace_dups <- function(dup_tip, origNodes, tipNamesOri, tip_names_toInsert, tip_names){

  #get index of duplicated IDs in tip_names vector
  index_to_remove <- which(tip_names %in% dup_tip)

  # save the values of tipNamesOri in a temporary object
  tip_tmp <- tipNamesOri

  tipNamesOri <- tip_tmp[-index_to_remove]
  # add new tip names
  # index to add tip name(s) in relation to the original vector of vertex IDs
  # origNodes does not have any duplicates
  index_to_add <- match(dup_tip, origNodes)
  # subtract 1 of the index to use with function append
  # append is an R function that will append a value after the index value of
  # interest
  index_to_add2 <- index_to_add - 1
  # This if will only check for consistency. If this does not hold, I coded
  # something not right.
  if(length(tip_names_toInsert) != length(index_to_add2)){
    stop("length of `tip_names_toInsert` should be the same as `index_to_add2`")
  }
  #sink( file = 'test.txt' , append = TRUE)
  while(length(index_to_add2) > 0){
    #print("tipNamesOri")
    #print(tipNamesOri[1:100])
    tipNamesOri2 <- append(x = tipNamesOri, values = tip_names_toInsert[1], after = index_to_add2[1])
    #print("tipNamesOri2")
    #print(tipNamesOri2[1:100])
    #print("origNodes")
    #print(origNodes[1:100])
    #print("tip_names_toInsert[1]")
    #print(tip_names_toInsert[1])
    #print("index_to_add2[1]")
    #print(index_to_add2[1])
    tipNamesOri <- tipNamesOri2
    index_to_add2 <- index_to_add2[-1]
    tip_names_toInsert <- tip_names_toInsert[-1]

  }
  #sink()

  return(tipNamesOri2)
}


#' Scale branch length from a phylogenetic tree
#'
#' @param tree a phylo object
#' @param scale numeric value to convert branch length
#'
#' @return a phylogenetic tree with the branch lengths
#'   scaled to vale of scale.
#' @export
convert_branches <- function(tree, scale = 1/365){

  #convert branch lengths from days to years ----
  tree$edge.length <- tree$edge.length * scale

  if(!is.null(tree$root.edge)){
    tree$root.edge <- tree$root.edge * scale
  }

  return(tree)
}


#' Reorder tip names
#'
#' Function to reorder the tip names of a phylogenetic tree.
#'
#' @param tip_names_migrant Vector of tip names in the form of ID_migrant
#' @param tip_names_vts Vector of tip names as returned by using the VirusTreeSimulator.
#'
#' @details Vector of tip_names_migrant will be reordered according to
#'    vector tip_names_vts.
#'
#' @return A vector with the reordered tip names
#' @export
reorder_tip_names <- function(tip_names_migrant, tip_names_vts){
  # Get tip names in the form of ID_migrant
  tip_names_migrant_ID <- str_match(string = tip_names_migrant, pattern = "\\d+")
  tip_names_vts_ID <- str_match(string = tip_names_vts, pattern = "\\d+")

  # reorder ID names
  # Saving indixes for how to reorder `tip_names_migrant_ID` to match `tip_names_vts_ID`
  reorder_idx <- match(tip_names_vts_ID, tip_names_migrant_ID)

  # Reordering the tip_names_migrant_ID  to match the order of the tip_names_ID vector
  reordered_tip_names <- tip_names_migrant[reorder_idx]

  return(reordered_tip_names)
}


#' Merge phylogenetic trees into a single tree
#'
#' Merge 2 or more phylogenetic trees into a single phylogenetic tree
#'
#' @param trees list of phylogenetic trees
#'
#' @return merged trees as object of class phylo
#' @export
merge_trees <- function(trees){


  #merge initial 2 trees
  treesxy <- bind.tree(trees[[1]], trees[[2]], position = trees[[1]]$root.edge)

  #number of trees to be added to the merged tree
  n_trees <- length(trees) - 2

  if(n_trees > 0){
    #index to merge next tree
    #index_count <- 3
    while(length(trees) - 2 > 0){
      tree <- trees[[3]]

      #add small branch length so as to correctly merge the next tree
      tree$root.edge <- tree$root.edge + 0.0001
      treesxy$root.edge <- 0.0001

      #merge tree
      new_tree <- bind.tree(treesxy, tree, position = treesxy$root.edge)
      treesxy <- new_tree
      trees <- trees[-3]
      #index_count <- index_count + 1
    }
  }

  return(treesxy)
}


#' Add root edge to phylogenetic tree
#'
#' @param tree object of class Phylo
#' @param total_sim_steps scalar for the total number of steps of network simulations
#' @param root.edge_value scalar for the value to use as root.edge in the phylogenetic tree
#'
#' @return A phylogenetic tree with a value for the root.edge.
#' @export
add_root_edge <- function(tree, total_sim_steps, root.edge_value = 7665){
  #browser()
  #max_edge <- max(distRoot(tree))
  max_edge <- max(get_all_distances_to_root(tree))
  tree_root.edge <- (total_sim_steps - max_edge) + root.edge_value
  tree$root.edge <- tree_root.edge

  return(tree)

}


#' Add root edge to phylogenetic tree
#'
#' @param tree object of class Phylo
#' @param root.edge_value scalar for the value to use as root.edge in the
#'    phylogenetic tree
#'
#' @return A phylogenetic tree with a value for root.edge.
#' @export
add_root_edge2 <- function(tree, root.edge_value = 0){
  #browser()
  #max_edge <- max(distRoot(tree))
  #max_edge <- max(get_all_distances_to_root(tree))
  #tree_root.edge <- (total_sim_steps - max_edge) + root.edge_value
  tree$root.edge <- root.edge_value

  return(tree)

}

#' Get tip names for cherries in a phylogenetic tree
#'
#' @param phy object of class phylo
#'
#' @return list of pairs of tip_names that are cherries
#' @export
get_tip_cherry <- function(phy){
  if (!inherits(phy, "phylo"))
    stop("object \"phy\" is not of class \"phylo\"")
  n <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != n - 1)
    stop("\"phy\" is not fully dichotomous")
  if (n < 4)
    stop("not enough tips in your phylogeny for this analysis")
  #get node number
  cherry_nodes <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n]) == 2)
  # get cherry node index
  cherry_index <- lapply(cherry_nodes, function(x) which(phy$edge[,1]==x))
  #get tip names
  tip_names <- lapply(cherry_index, function(x) phy$tip.label[phy$edge[x,2]])

  return(tip_names)
}
