% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling_utils.R
\name{sampleIDs2}
\alias{sampleIDs2}
\title{Sample IDs of individuals who have been diagnosed, are active and on or not on ART}
\usage{
sampleIDs2(
  perc,
  start_date,
  end_date,
  art_init,
  departure,
  diag_info,
  origin,
  tm,
  location
)
}
\arguments{
\item{perc}{Percentage of the population to be sampled: Numbers should be between
0 and 1.}

\item{start_date}{Start date for sampling times in decimal year.}

\item{end_date}{End date for sampling times in decimal year.}

\item{art_init}{Dataframe of IDs and time each ID initiated antiretroviral
treatment (ART). Time must be in decimal year.}

\item{departure}{Dataframe of IDs and their time of departure. Time must be
in decimal year.}

\item{diag_info}{Dataframe of IDs and their time of diagnosis. Time must be
in decimal year.}

\item{tm}{Transmission matrix of who transmitted to whom.}

\item{location}{Character of location at time of transmission to get IDs
to sample from. Location can take the value of "region" or "global".}
}
\value{
Dataframe of IDs and their sampled times
}
\description{
Function to sample IDs based on sampled time and metadata. This function will
sample individuals that are diagnosed and active (individuals can or cannot be on ART)
at the time of sampling.
}
\details{
This function will get all IDs in transmission matrix (tm) from
   location. From those IDs, it selected all IDs that have been diagnosed
   during a simulation run. From all diagnosed IDs, it will select all IDs
   that have been diagnosed before sampled time. Then, from this new list of
   IDs, it select all IDs that have been active before time of sampling.
   Finally, from active and diagnosed IDs before sampled times, it will get
   all IDs that never initiated ART before sample times. From those IDs not
   on ART, the function will random sample a ID to be associated to the
   sampled time. This function is to mimic the idea of sampling individuals
   on HIV studies.
}
