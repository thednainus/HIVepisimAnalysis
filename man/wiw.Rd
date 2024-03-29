% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_utils.R
\name{wiw}
\alias{wiw}
\title{Who infected whom}
\usage{
wiw(df1)
}
\arguments{
\item{df1}{dataframe containing information of the pairs that a TP (true positives),
and the correct direction of transmission.}
}
\value{
dataframe containg the absolute number of the total pairs that are
   TP and from those the total number of pairs in which phyloscanner correctly
   estimated who infected whom.
}
\description{
Who infected whom
}
\details{
This function will return the absolute numbers for the total number
   of pairs that are TP and from those pairs the total number in which
   phyloscanner also correctly estimated who infected whom.
}
