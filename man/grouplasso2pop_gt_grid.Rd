% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grouplasso2pop.R
\name{grouplasso2pop_gt_grid}
\alias{grouplasso2pop_gt_grid}
\title{Compute group lasso for two populations with group testing data over a grid of tuning parameter values}
\usage{
grouplasso2pop_gt_grid(
  Y1,
  Z1,
  Se1,
  Sp1,
  X1,
  groups1,
  E.approx1 = FALSE,
  Y2,
  Z2,
  Se2,
  Sp2,
  X2,
  groups2,
  E.approx2 = FALSE,
  rho1,
  rho2,
  n.lambda,
  n.eta,
  lambda.min.ratio,
  lambda.max.ratio = 1,
  eta.min.ratio = 0.001,
  eta.max.ratio = 10,
  w1,
  w2,
  w,
  AA1,
  AA2,
  Com,
  tol = 0.001,
  maxiter = 1000,
  report.prog = TRUE
)
}
\arguments{
\item{Y1}{Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Z1}{Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Se1}{A vector of testing sensitivities, where the first element is the
testing specificity for pools and the second entry is the
test specificity for individual testing, if applicable.}

\item{Sp1}{A vector of testing specificities, where the first element is the
testing specificity for pools and the second entry is the
test specificity for individual testing, if applicable.}

\item{X1}{the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)}

\item{groups1}{a vector indicating to which group each covariate of data set 2 belongs}

\item{E.approx1}{a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly for data set 1}

\item{Y2}{Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Z2}{Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Se2}{A vector of testing sensitivities, where the first element is the
testing specificity for pools and the second entry is the
test specificity for individual testing, if applicable.}

\item{Sp2}{A vector of testing specificities, where the first element is the
testing specificity for pools and the second entry is the
test specificity for individual testing, if applicable.}

\item{X2}{the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)}

\item{groups2}{a vector indicating to which group each covariate of data set 2 belongs}

\item{E.approx2}{a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly for data set 2}

\item{rho1}{weight placed on the first data set}

\item{rho2}{weight placed on the second data set}

\item{n.lambda}{the number of lambda values}

\item{n.eta}{the number of eta values}

\item{lambda.min.ratio}{ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model}

\item{lambda.max.ratio}{ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model}

\item{eta.min.ratio}{ratio of the smallest to largest value in the sequence of eta values}

\item{eta.max.ratio}{controls the largest value of eta in the eta sequence}

\item{w1}{group-specific weights for different penalization across groups in data set 1}

\item{w2}{group-specific weights for different penalization across groups in data set 2}

\item{w}{group-specific weights for different penalization toward similarity for different groups}

\item{AA1}{a list of the matrices A2j}

\item{Com}{the indices of the covariate groups which are common in the two datasets}

\item{tol}{a convergence criterion}

\item{maxiter}{the maximum allowed number of iterations (EM steps)}

\item{report.prog}{a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.}
}
\value{
Returns the estimator of the parametric model with group testing data
}
\description{
Compute group lasso for two populations with group testing data over a grid of tuning parameter values
}
\examples{
grouplasso2pop_gt_data <- get_grouplasso2pop_data( n1 = 1000, n2 = 1200, response = "gt")

grouplasso2pop_gt_grid.out <- grouplasso2pop_gt_grid(Y1 = grouplasso2pop_gt_data$Y1$I,
                                                     Z1 = grouplasso2pop_gt_data$Y1$A,
                                                     Se1 = grouplasso2pop_gt_data$Y1$Se,
                                                     Sp1 = grouplasso2pop_gt_data$Y1$Sp,
                                                     X1 = grouplasso2pop_gt_data$X1,
                                                     groups1 = grouplasso2pop_gt_data$groups1,
                                                     E.approx1 = grouplasso2pop_gt_data$Y1$E.approx,
                                                     Y2 = grouplasso2pop_gt_data$Y2$I,
                                                     Z2 = grouplasso2pop_gt_data$Y2$A,
                                                     Se2 = grouplasso2pop_gt_data$Y2$Se,
                                                     Sp2 = grouplasso2pop_gt_data$Y2$Sp,
                                                     X2 = grouplasso2pop_gt_data$X2,
                                                     groups2  = grouplasso2pop_gt_data$groups2,
                                                     E.approx2 = grouplasso2pop_gt_data$Y2$E.approx,
                                                     rho1 = 1,
                                                     rho2 = 1,
                                                     n.lambda = 10,
                                                     n.eta = 5,
                                                     lambda.min.ratio = 0.01,
                                                     lambda.max.ratio = 0.50,
                                                     eta.min.ratio = 0.01,
                                                     eta.max.ratio = 1,
                                                     w1 = grouplasso2pop_gt_data$w1,
                                                     w2 = grouplasso2pop_gt_data$w2,
                                                     w = grouplasso2pop_gt_data$w,
                                                     AA1 = grouplasso2pop_gt_data$AA1,
                                                     AA2 = grouplasso2pop_gt_data$AA2,
                                                     Com = grouplasso2pop_gt_data$Com,
                                                     tol = 1e-3,
                                                     maxiter = 500,
                                                     report.prog = TRUE)
}
