% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_utils.R
\name{get_rates}
\alias{get_rates}
\title{Get rates to construct ROC curve}
\usage{
get_rates(threshold, df_true)
}
\arguments{
\item{threshold}{list object with threshold values rangng from 0 to 1}

\item{df_true}{data.frame object for the infector probability and the true
classification of transmission pairs}
}
\value{
data.frame object for values of threshold, true positive rate (TPR)
   and false positive rate (FPR).
}
\description{
Get false positive rates and true positive rate for ROC curves
}
\details{
The true data for df_true object will have 1 for a transmission pair
   that have occurred and 0 for a transmission pair that did not occur. This
   will be based by comparison with the transmission matrix.
}
