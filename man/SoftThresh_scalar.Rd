% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SoftThresh_scalar}
\alias{SoftThresh_scalar}
\title{The soft-thresholding function for a scalar input}
\usage{
SoftThresh_scalar(z, a)
}
\arguments{
\item{z}{the value to which to apply soft-thresholding}

\item{a}{the threshold}
}
\value{
the value of the soft-thresholding function
}
\description{
The soft-thresholding function for a scalar input
}
\examples{
z <- 3
a <- 2
SoftThresh(z,a)
}
