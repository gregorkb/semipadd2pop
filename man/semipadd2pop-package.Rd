\name{semipadd2pop-package}
\alias{semipadd2pop-package}
\alias{semipadd2pop}
\docType{package}
\title{
  Combine two data sets to fit semiparametric additive models
}
\description{
\packageDescription{semipadd2pop}
  
}
\details{
  The principal functions in the package for fitting sparse semiparametric regression models on single data set are \code{semipadd} and \code{semipadd_cv_adapt}.  The first requires the user to choose the sparsity tuning parameter value, while the second chooses it via crossvalidation and includes an adaptive step. The functions \code{semipadd2pop} and \code{semipadd2pop_cv_adapt} are for fitting sparse semiparametric regression models on two data sets which have some covariates in common, such that the fitted effects of the common covariates are encouraged via a ridge-like penalty to be similar.  The first requires the user to choose all tuning parameters, while the second selects the tuning parameters governing sparsity and penalization towards similarity via crossvalidation, using an adaptive step.
}
\author{
\packageAuthor{semipadd2pop}

Maintainer: \packageMaintainer{semipadd2pop}
}
\references{
  This is a package for a working paper with Chris McMahan and Dewei Wang entitled ``Combining data sets for semiparametric additive modeling''.
}
\keyword{ package }
\examples{
  \dontrun{
# fit a semiparametric model on a single simulated data set
data <- get_semipadd_data(n = 500, response = "continuous")

semipadd.out <- semipadd(Y = data$Y,
                         X = data$X,
                         nonparm = data$nonparm,
                         response = "continuous",
                          w = 1,
                          d = 20,
                          xi = 1,
                          lambda.beta = 1,
                          lambda.f = 1,
                          tol = 1e-3,
                          max.iter = 500)
 
plot_semipadd(semipadd.out, 
              true.functions = list( f = data$f,
                                     X = data$X))

# fit semiparametric models combining two simulated data sets with some common covariates                                    
data <- get_semipadd2pop_data(n1 = 500, n2 = 600, response = "binary")

semipadd2pop.out <- semipadd2pop(Y1 = data$Y1,
                                 X1 = data$X1,
                                 nonparm1 = data$nonparm1,
                                 Y2 = data$Y2,
                                 X2 = data$X2,
                                 nonparm2 = data$nonparm2,
                                 response = "binary",
                                 rho1 = 1,
                                 rho2 = 1,
                                 nCom = data$nCom,
                                 d1 = 25,
                                 d2 = 15,
                                 xi = .5,
                                 w1 = 1,
                                 w2 = 1,
                                 w = 1,
                                 lambda.beta = .01,
                                 lambda.f = .01,
                                 eta.beta = .1,
                                 eta.f = .1,
                                 tol = 1e-3,
                                 maxiter = 500)
                             
plot_semipadd2pop(semipadd2pop.out,
                  true.functions = list(f1 = data$f1,
                                        f2 = data$f2,
                                        X1 = data$X1,
                                        X2 = data$X2))
               
  }
}
