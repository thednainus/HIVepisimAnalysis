% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_utils.R
\name{get_chain}
\alias{get_chain}
\title{Get transmission chains}
\usage{
get_chain(df1, tm)
}
\arguments{
\item{df1}{Dataframe containing information of summarized results obtained
with phyloscanner after running the function summarize_trans.}

\item{tm}{Transmission matrix containing the information of who transmitted
whom.}

\item{all_trans}{List of dataframes containing the transmission chains for
each ID.}
}
\value{
A dataframe with the transmission chain and the time of transmissions
}
\description{
Get transmission chains
}
\details{
This function will go over each transmission to create a "transmission
chain from ID1 to ID5, for example. This will help indentify in which
circunstances phyloscanner indetified two IDs as a transmission pair,
when they were not one.
}
