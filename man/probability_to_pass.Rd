% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equating.R
\name{probability_to_pass}
\alias{probability_to_pass}
\title{The probability to pass on a reference test given a score on a new booklet}
\usage{
probability_to_pass(
  dataSrc,
  parms,
  ref_items,
  pass_fail,
  predicate = NULL,
  target_booklets = NULL,
  nDraws = 1000
)
}
\arguments{
\item{dataSrc}{a connection to a dexter database, a matrix, or a data.frame with columns: person_id, item_id, item_score}

\item{parms}{object produced by \code{\link{fit_enorm}} or a data.frame with columns item_id, item_score and  beta. 
If uncertainty about parameter estimation should be included
in the computations, use a `parms` object computed with `method='Bayes'` and nDraws equal or larger than nDraws in probability_to_pass}

\item{ref_items}{vector with id's of items in the reference set, they must all occur in dataSrc}

\item{pass_fail}{pass-fail score on the reference set, the lowest score with which one passes}

\item{predicate}{An optional expression to subset data in dataSrc, if NULL all data is used}

\item{target_booklets}{The target test booklet(s). A data.frame with columns booklet_id (if multiple booklets) and item_id, 
if NULL (default) this will be derived from the dataSrc and the probability to pass will be computed 
for each test score for each booklet in your data.}

\item{nDraws}{The function uses an Markov-Chain Monte-Carlo method to calculate the probability to pass and this is the number of Monte-Carlo samples used.}
}
\value{
An object of type \code{p2pass}. Use \code{coef()} to extract the 
probablity to pass for each booklet and score. Use \code{plot()} to plot 
the probabilities, sensitivity and specificity or a ROC-curve.
}
\description{
Given response data that form a connected design,
compute the probability to pass on the reference set conditional on each score on one or more target tests.
}
\details{
Note that this function is computationally intensive and can take some time to run, especially when computing the
probability to pass for multiple target booklets. Further technical details can be found in a vignette.
}
\seealso{
The function used to plot the results: \code{\link{plot.p2pass}}
}
