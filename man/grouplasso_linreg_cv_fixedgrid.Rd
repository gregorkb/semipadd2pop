% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grouplasso.R
\name{grouplasso_linreg_cv_fixedgrid}
\alias{grouplasso_linreg_cv_fixedgrid}
\title{Choose tuning parameters by crossvalidation for grouplasso linreg when given a fixed grid of lambda values}
\usage{
grouplasso_linreg_cv_fixedgrid(
  Y,
  X,
  groups,
  lambda.seq,
  n.folds,
  b.init.mat,
  w,
  tol = 0.001,
  maxiter = 500
)
}
\arguments{
\item{Y}{the response vector}

\item{X}{matrix containing the design matrices}

\item{groups}{a vector indicating to which group each covariate belongs}

\item{lambda.seq}{sequence of lambda values}

\item{n.folds}{the number of crossvalidation folds}

\item{w}{group-specific weights for different penalization of different groups}

\item{tol}{the convergence tolerance}

\item{maxiter}{the maximum number of iterations allowed for each fit}

\item{b.init.arr}{array of initial values for beta}
}
\value{
a list containing the fits over a grid of lambda values as well as the vector of lambda values
}
\description{
Choose tuning parameters by crossvalidation for grouplasso linreg when given a fixed grid of lambda values
}
\examples{
grouplasso_linreg_data <- get_grouplasso_data(n = 400, response = "continuous")

grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = grouplasso_linreg_data$Y,
                                                     X = grouplasso_linreg_data$X,
                                                     groups = grouplasso_linreg_data$groups,
                                                     n.lambda = 25,
                                                     lambda.min.ratio = 0.001,
                                                     lambda.max.ratio = 0.1,
                                                     w = grouplasso_linreg_data$w,
                                                     tol = 1e-3,
                                                     maxiter = 500,
                                                     report.prog = FALSE)

grouplasso_linreg_cv_fixedgrid.out <- grouplasso_linreg_cv_fixedgrid(Y = grouplasso_linreg_data$Y,
                                                                     X = grouplasso_linreg_data$X,
                                                                     groups = grouplasso_linreg_data$groups,
                                                                     lambda.seq = grouplasso_linreg_grid.out$lambda.seq,
                                                                     n.folds = 5,
                                                                     b.init.mat = grouplasso_linreg_grid.out$b.mat,
                                                                     w = grouplasso_linreg_data$w,
                                                                     tol = 1e-3,
                                                                     maxiter = 500)
}
