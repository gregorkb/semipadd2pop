% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{EYapprox}
\alias{EYapprox}
\title{Approximates the conditional expectations of individual disease statuses with Gibbs sampling.}
\usage{
EYapprox(Z, Y, X, b, Se, Sp, GI = 5000)
}
\arguments{
\item{Z}{Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Y}{Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{X}{Design matrix with first column a column of 1s.}

\item{b}{Parameter values at which to compute the conditional expectations.}

\item{Se}{A vector of testing sensitivities of length \code{max(Z[,3])}.}

\item{Sp}{A vector of testing specificities of length \code{max(Z[,3])}.}

\item{GI}{The length of the Gibbs sampling Markov chain.}
}
\value{
The vector of conditional expectations.

This function uses a Gibbs sampler to appriximate the conditional expectation of 
each individual's disease status, conditional on the observed assay data and the disease 
statuses of all other individuals. This function is used in the EM algorithm
performed by the functions \code{mlegt}, \code{enetgt}, \code{enetgt.grid}, and 
\code{enetgt.grid.0} under array testing.
}
\description{
Approximates the conditional expectations of individual disease statuses with Gibbs sampling.
}
\examples{
# generate individual covariate values and disease statuses
grouplasso_gt_data <- get_grouplasso_gt_data(100)
EY <- EYapprox(Z = grouplasso_gt_data$Z,
               Y = grouplasso_gt_data$Y,
               X = grouplasso_gt_data$X,
               b = grouplasso_gt_data$beta,
               Se = grouplasso_gt_data$Se,
               Sp = grouplasso_gt_data$Sp,
               GI = 10000)
}
