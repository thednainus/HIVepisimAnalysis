% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/W_utils.R
\name{validadeW}
\alias{validadeW}
\title{Title Validate infector probability calculation}
\usage{
validadeW(
  sim,
  run,
  tm,
  W_true,
  W_estimated,
  seed_ID,
  MH,
  true_tree,
  prefix = NULL
)
}
\arguments{
\item{sim}{Simulation number to know the parameter values that network was
simulated}

\item{run}{Run number}

\item{tm}{Transmission matrix per seed}

\item{W_true}{Infector probability calculated on the true trees}

\item{seed_ID}{Seed ID name. Seed is the node ID that started the
transmissions; it is the initial infected node in the simulations.}

\item{MH}{maximum height in which infector probabilities were calculated}

\item{true_tree}{object of class phylo}

\item{prefix}{prefix to save results. If prefix = NULL, results will be saved
on file W_stats.csv.}
}
\value{
a data frame of 1 row with the following values:
 1. Code = code indicating if W was calculated on true trees
 2. Maximum height = maximum height (MH) in which calculations were performed;
 3. Total tips = total number of tips in the phylogenetic tree;
 4. n_tips_region = total number of tips from region in the phylogenetic tree;
 5. n_tips_global = total number of tips from global in the phylogenetic tree;
 6. n_trans_W = total number of transmissions independent of infector probability
    values;
 7. n_trans_tm = total number of transmissions based on transmission matrix
    subset by MH;
 8. n_true_trans_W = total number of true transmission independent of infector
    probability (W) values. Here is the true transmission when compared to the
    transmission matrix;
 9. n_W80_all = number of transmission in which W >= 80%;
 10. n_W80_correctDonorRecipt = number of transmission in point 9 that is a
     true transmission (when comparing to the transmission matrix);
 11. n_W80_swapDonorRecipt = number of transmission in point 9 that is a
     true transmission (when comparing to the transmission matrix), but
     there was a swap between donor and recipient;
 12. expected_n_trans = expected number of transmission defined as the sum
     of all infector probability (W) values.
}
\description{
This function saves several values to compare W calculated on the true trees
   and the true transmission matrix. For the return values see below.
}
