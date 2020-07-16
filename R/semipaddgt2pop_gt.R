
#' Compute group lasso for two populations with group testing data
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda the level of sparsity penalization
#' @param eta the level of penalization towards model similarity
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param init a list of initial values for the coefficient for each of the two data sets
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#' @examples 
#' grouplassogt2pop_data <- get_grouplassogt2pop_data( n1 = 1000, n2 = 1200)
#' 
#' grouplassogr2pop.out <- grouplassogt2pop(Y1 = grouplassogt2pop_data$Y1,
#'                                          Z1 = grouplassogt2pop_data$Z1,
#'                                          Se1 = grouplassogt2pop_data$Se1,
#'                                          Sp1 = grouplassogt2pop_data$Sp1,
#'                                          X1 = grouplassogt2pop_data$X1,
#'                                          groups1 = grouplassogt2pop_data$groups1,
#'                                          Y2 = grouplassogt2pop_data$Y2,
#'                                          Z2 = grouplassogt2pop_data$Z2,
#'                                          Se2 = grouplassogt2pop_data$Se2,
#'                                          Sp2 = grouplassogt2pop_data$Sp2,
#'                                          X2 = grouplassogt2pop_data$X2,
#'                                          groups2  = grouplassogt2pop_data$groups2,
#'                                          rho1 = 2,
#'                                          rho2 = 1,
#'                                          lambda = 1,
#'                                          eta = 1,
#'                                          w1  = grouplassogt2pop_data$w1,
#'                                          w2  = grouplassogt2pop_data$w2,
#'                                          w  = grouplassogt2pop_data$w,
#'                                          AA1  = grouplassogt2pop_data$AA1,
#'                                          AA2  = grouplassogt2pop_data$AA2,
#'                                          Com  = grouplassogt2pop_data$Com,
#'                                          tol = 1e-3,
#'                                          max.iter = 500,
#'                                          report.prog = TRUE)
#' @export
grouplassogt2pop <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,lambda,eta,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,init=NULL,report.prog=TRUE)
{

  # set function to get conditional expectations
  get.EY <- eval(parse(text = ifelse(E.approx, "EYapprox","EYexact")))
  
  # set initial values
  if(length(init) == 0){

    beta1.hat1 <- rep(0,ncol(X1))
    beta2.hat1 <- rep(0,ncol(X2))

  } else{

    beta1.hat1 <- init$beta1
    beta2.hat1 <- init$beta2

  }

  ###### Do the EM-algorithm with penalized updates
  conv <- 1
  iter <- 0
  inner.iter <- numeric()
  while( conv > tol & iter < max.iter)
  {

    beta1.hat0 <- beta1.hat1
    beta2.hat0 <- beta2.hat1

    # E-step: compute the conditional expectations for the true disease statuses
    EY1 <- as.numeric(get.EY(Z1,Y1,X=X1,b=beta1.hat1,Se1,Sp1))
    EY2 <- as.numeric(get.EY(Z2,Y2,X=X2,b=beta2.hat1,Se2,Sp2))

    # update initial values
    init <- list(beta1 = beta1.hat1,
                 beta2 = beta2.hat1)

    # M-step: maximize the objective function with conditional expectations substituted
    grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = EY1,
                                                       rX1 = X1,
                                                       groups1 = groups1,
                                                       rY2 = EY2,
                                                       rX2 = X2,
                                                       groups2 = groups2,
                                                       rho1 = rho1,
                                                       rho2 = rho2,
                                                       lambda = lambda,
                                                       eta = eta,
                                                       w1 = w1,
                                                       w2 = w2,
                                                       w = w,
                                                       rAA1 = AA1,
                                                       rAA2 = AA2,
                                                       rCom = Com,
                                                       tol = tol,
                                                       maxiter = max.iter,
                                                       beta1_init = init$beta1,
                                                       beta2_init = init$beta2)

    beta1.hat1 <- grouplasso2pop_logreg.out$beta1.hat
    beta2.hat1 <- grouplasso2pop_logreg.out$beta2.hat

    conv <- max(abs(c(beta1.hat1 - beta1.hat0, beta2.hat1 - beta2.hat0)))
    iter <- iter + 1
    inner.iter[iter] <- grouplasso2pop_logreg.out$iter

    if(report.prog) print(grouplasso2pop_logreg.out$iter)

  }

  output <- list( beta1.hat = beta1.hat1,
                  beta2.hat = beta2.hat1,
                  inner.iter = inner.iter)

}

#' Compute group lasso for two populations with group testing data over a grid of tuning parameter values
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values
#' @param n.eta the number of eta values
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#' @export
grouplassogt2pop_grid <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE)
{
  # pull the individual diagnoses to be treated as true disease statuses to find the sequence of lambda values.
  # later on think about how to do this without individual diagnoses, e.g. under master pool testing.
  Y1.diag <- pull.diagnoses(Z1,Y1)
  Y2.diag <- pull.diagnoses(Z2,Y2)

  # find lambda.max
  q1 <- length(unique(groups1))
  q2 <- length(unique(groups2))

  n1 <- nrow(X1)
  n2 <- nrow(X2)

  norms1 <- numeric(q1)
  norms2 <- numeric(q2)
  for(j in 2:max(q1,q2))
  {

    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)

    if(j <= q1){
      
      norms1[j] <- sqrt(sum((t(X1[,ind1]) %*% (Y1.diag - mean(Y1.diag)))^2)) / ( w1[j] * n1 / rho1)
      
    }
    if(j <= q2){
      
      norms2[j] <- sqrt(sum((t(X2[,ind2]) %*% (Y2.diag - mean(Y2.diag)))^2)) / (w2[j] * n2 / rho1)
      
    }

  }

  lambda.max <- 2 * max(norms1,norms2)

  # make a lambda sequence
  lambda.min <- lambda.min.ratio * lambda.max
  lambda.seq <- sort(c(exp(log(lambda.min) + ((n.lambda+1):1)/(n.lambda+1) * ((log(lambda.max) - log(lambda.min)))))[-1])

  if(n.lambda == 1) lambda.seq <- lambda.min

  # make the eta sequence after fitting the model for the smallest value of lambda
  eta.seq <- numeric(n.eta)

  b1.arr <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  b2.arr <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  iterations <- matrix(0,n.lambda*n.eta,3)
  colnames(iterations) <- c("lambda","eta","EM-steps")
  steps <- 0
  init <- NULL
  for(l in 1:n.lambda){
    for(k in 1:n.eta){

      grouplassogt2pop.out <- grouplassogt2pop(Y1 = Y1,
                                               Z1 = Z1,
                                               Se1 = Se1,
                                               Sp1 = Sp1,
                                               X1 = X1,
                                               groups1 = groups1,
                                               Y2 = Y2,
                                               Z2 = Z2,
                                               Se2 = Se2,
                                               Sp2 = Sp2,
                                               X2 = X2,
                                               groups2 = groups2,
                                               rho1 = rho1,
                                               rho2 = rho2,
                                               lambda = lambda.seq[l],
                                               eta = eta.seq[k],
                                               w1 = w1,
                                               w2 = w2,
                                               w = w,
                                               AA1 = AA1,
                                               AA2 = AA2,
                                               Com = Com,
                                               E.approx = E.approx,
                                               tol = tol,
                                               max.iter = max.iter,
                                               init = init,
                                               report.prog=FALSE)

      b1 <- grouplassogt2pop.out$beta1.hat
      b2 <- grouplassogt2pop.out$beta2.hat

      if(l == 1 & k == 1){# define the eta sequence

        P1 <- logit(X1 %*% b1)
        P2 <- logit(X2 %*% b2)

        neg2LL1 <- - 2 * rho1 / n1 * sum(Y1.diag*log(P1) + (1 - Y1.diag)*log(1-P1))
        neg2LL2 <- - 2 * rho2 / n2 * sum(Y2.diag*log(P2) + (1 - Y2.diag)*log(1-P2))

        beta1beta2.wl2 <- 0
        for(j in Com)
        {
          ind1 <- which(groups1 == j)
          ind2 <- which(groups2 == j)
          beta1beta2.wl2  <- beta1beta2.wl2 + w[j] * sum( (AA1[[j]] %*% b1[ind1] - AA2[[j]] %*% b2[ind2] )^2 )
        }

        eta.max <- (neg2LL1 + neg2LL2) / beta1beta2.wl2
        eta.min <- 0.001 * eta.max
        eta.seq <- sort(exp(log(eta.min) + ((n.eta-1):0)/(n.eta-1) * (log(eta.max) - log(eta.min))) )
        eta.seq[1] <- 0

      }

      init <- list( beta1 = b1,
                    beta2 = b2)

      b1.arr[,l,k] <- b1
      b2.arr[,l,k] <- b2

      steps <- steps + 1
      iterations[steps,] <- c(lambda.seq[l],eta.seq[k],length(grouplassogt2pop.out$inner.iter))

      if(report.prog == TRUE){

        print(c(l,k,length(grouplassogt2pop.out$inner.iter)))

      }

    }

  }

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = iterations)
}


#' Compute group lasso for two populations with group testing data over a user-specified grid of tuning parameter values
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda.seq the sequence of lambda values
#' @param eta.seq the sequence of eta values
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#' @export
grouplassogt2pop_fixedgrid <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE)
{

  n.lambda <- length(lambda.seq)
  n.eta <- length(eta.seq)

  b1.arr <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  b2.arr <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  iterations <- matrix(0,n.lambda*n.eta,3)
  colnames(iterations) <- c("lambda","eta","EM-steps")
  steps <- 0
  init <- NULL
  for(l in 1:n.lambda){
    for(k in 1:n.eta){

      grouplassogt2pop.out <- grouplassogt2pop(Y1 = Y1,
                                               Z1 = Z1,
                                               Se1 = Se1,
                                               Sp1 = Sp1,
                                               X1 = X1,
                                               groups1 = groups1,
                                               Y2 = Y2,
                                               Z2 = Z2,
                                               Se2 = Se2,
                                               Sp2 = Sp2,
                                               X2 = X2,
                                               groups2 = groups2,
                                               rho1 = rho1,
                                               rho2 = rho2,
                                               lambda = lambda.seq[l],
                                               eta = eta.seq[k],
                                               w1 = w1,
                                               w2 = w2,
                                               w = w,
                                               AA1 = AA1,
                                               AA2 = AA2,
                                               Com = Com,
                                               E.approx = E.approx,
                                               tol = tol,
                                               max.iter = max.iter,
                                               init = init,
                                               report.prog = FALSE)

      b1 <- grouplassogt2pop.out$beta1.hat
      b2 <- grouplassogt2pop.out$beta2.hat

      init <- list( beta1 = b1,
                    beta2 = b2)

      b1.arr[,l,k] <- b1
      b2.arr[,l,k] <- b2

      steps <- steps + 1
      iterations[steps,] <- c(lambda.seq[l],eta.seq[k],length(grouplassogt2pop.out$inner.iter))

      if(report.prog == TRUE){

        print(c(l,k,length(grouplassogt2pop.out$inner.iter)))

      }

    }

  }

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = iterations)
}
#' Choose tuning parameters by crossvalidation for grouplassogt2pop when given a fixed grid of lambda and eta values
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda.seq the sequence of lambda values
#' @param eta.seq the sequence of eta values
#' @param n.folds the number of crossvalidation folds
#' @param b1.init.arr array of initial values for beta1
#' @param b2.init.arr array of initial values for beta2
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#'
#' @examples
#' grouplassogt2pop_data <- get_grouplassogt2pop_data(n1=400,n2=600)
#'
#' grouplassogt2pop_grid.out <- grouplassogt2pop_grid(Y1 = grouplassogt2pop_data$Y1,
#'                                                    Z1 = grouplassogt2pop_data$Z1,
#'                                                    Se1 = grouplassogt2pop_data$Se1,
#'                                                    Sp1 = grouplassogt2pop_data$Sp1,
#'                                                    X1 = grouplassogt2pop_data$X1,
#'                                                    groups1 = grouplassogt2pop_data$groups1,
#'                                                    Y2 = grouplassogt2pop_data$Y2,
#'                                                    Z2 = grouplassogt2pop_data$Z2,
#'                                                    Se2 = grouplassogt2pop_data$Se2,
#'                                                    Sp2 = grouplassogt2pop_data$Sp2,
#'                                                    X2 = grouplassogt2pop_data$X2,
#'                                                    groups2 = grouplassogt2pop_data$groups2,
#'                                                    rho1 = 2,
#'                                                    rho1 = 1,
#'                                                    n.lambda = 5,
#'                                                    n.eta = 5,
#'                                                    lambda.min.ratio = 0.01,
#'                                                    w1 = grouplassogt2pop_data$w1,
#'                                                    w2 = grouplassogt2pop_data$w2,
#'                                                    w = grouplassogt2pop_data$w,
#'                                                    AA1 = grouplassogt2pop_data$AA1,
#'                                                    AA2 = grouplassogt2pop_data$AA2,
#'                                                    Com = grouplassogt2pop_data$Com,
#'                                                    tol = 1e-2,
#'                                                    max.iter = 500,
#'                                                    report.prog = TRUE)
#'
#' grouplassogt2pop_cv_fixedgrid.out <- grouplassogt2pop_cv_fixedgrid(Y1 = grouplassogt2pop_data$Y1,
#'                                                                    Z1 = grouplassogt2pop_data$Z1,
#'                                                                    Se1 = grouplassogt2pop_data$Se1,
#'                                                                    Sp1 = grouplassogt2pop_data$Sp1,
#'                                                                    X1 = grouplassogt2pop_data$X1,
#'                                                                    groups1 = grouplassogt2pop_data$groups1,
#'                                                                    Y2 = grouplassogt2pop_data$Y2,
#'                                                                    Z2 = grouplassogt2pop_data$Z2,
#'                                                                    Se2 = grouplassogt2pop_data$Se2,
#'                                                                    Sp2 = grouplassogt2pop_data$Sp2,
#'                                                                    X2 = grouplassogt2pop_data$X2,
#'                                                                    groups2 = grouplassogt2pop_data$groups2,
#'                                                                    rho1 = 2,
#'                                                                    rho2 = 1,
#'                                                                    lambda.seq = grouplassogt2pop_grid.out$lambda.seq,
#'                                                                    eta.seq = grouplassogt2pop_grid.out$eta.seq,
#'                                                                    n.folds = 5,
#'                                                                    b1.init.arr = grouplassogt2pop_grid.out$b1.arr,
#'                                                                    b2.init.arr = grouplassogt2pop_grid.out$b2.arr,
#'                                                                    w1 = grouplassogt2pop_data$w1,
#'                                                                    w2 = grouplassogt2pop_data$w2,
#'                                                                    w = grouplassogt2pop_data$w,
#'                                                                    AA1 = grouplassogt2pop_data$AA1,
#'                                                                    AA2 = grouplassogt2pop_data$AA2,
#'                                                                    Com = grouplassogt2pop_data$Com,
#'                                                                    tol = 1e-2,
#'                                                                    max.iter = 500)
#' @export
grouplassogt2pop_cv_fixedgrid <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,n.folds,b1.init.arr,b2.init.arr,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=500)
{

  # for now just use the diagnoses as the true responses and use the logreg fits on these
  Y1.diag <- pull.diagnoses(Z1,Y1)
  Y2.diag <- pull.diagnoses(Z2,Y2)

  # create list of sets of indices indicating which observations are in each fold
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  folds1 <- vector("list",n.folds)
  folds2 <- vector("list",n.folds)
  fold.size1 <- floor(n1 / n.folds)
  fold.size2 <- floor(n2 / n.folds)
  for(fold in 1:n.folds){

    folds1[[fold]] <- ((fold-1)*fold.size1 + 1):(fold*fold.size1)
    folds2[[fold]] <- ((fold-1)*fold.size2 + 1):(fold*fold.size2)
  }

  if( floor(n1 / n.folds) != n1/n.folds )
  {
    folds1[[n.folds]] <- c(folds1[[n.folds]],(fold*fold.size1+1):n1)
  }

  if( floor(n2 / n.folds) != n2/n.folds )
  {
    folds2[[n.folds]] <- c(folds2[[n.folds]],(fold*fold.size2+1):n2)
  }
  # get fits at all lambda and eta combinations on all cv folds

  n.lambda <- length(lambda.seq)
  n.eta <- length(eta.seq)

  b1.folds.arr <- array(0,dim=c(ncol(X1),n.lambda,n.eta,n.folds))
  b2.folds.arr <- array(0,dim=c(ncol(X2),n.lambda,n.eta,n.folds))

  minus2ll.arr <- array(0,dim=c(n.lambda,n.eta,n.folds))

  iterations <- matrix(0,n.lambda*n.eta,2+n.folds)
  colnames(iterations) <- c("lambda","eta",paste("fold",1:n.folds,"iter"))
  step <- 1

  for(l in 1:n.lambda){
    for(k in 1:n.eta){

      iterations[step,c(1,2)] <- c(lambda.seq[l],eta.seq[k])

      for(fold in 1:n.folds){

        fold1 <- folds1[[fold]]
        fold2 <- folds2[[fold]]

        # later on consider a "real" crossvalidation, where the folds consist of some master pools, as in aenetgt
        # grouplassogt2pop.out <- grouplassogt2pop(Y1,
        #                                          Z1,
        #                                          Se1 = Se1,
        #                                          Sp1 = Sp1,
        #                                          X1,
        #                                          groups1 = groups1,
        #                                          Y2,
        #                                          Z2,
        #                                          Se2 = Se2,
        #                                          Sp2 = Sp2,
        #                                          X2,
        #                                          groups2 = groups2,
        #                                          lambda = lambda.seq[l]*(n.folds - 1)/n.folds,
        #                                          eta = eta.seq[k]*(n.folds - 1)/n.folds,
        #                                          w1 = w1,
        #                                          w2 = w2,
        #                                          w = w,
        #                                          AA1 = AA1,
        #                                          AA2 = AA2,
        #                                          Com = Com,
        #                                          E.approx = E.approx,
        #                                          tol = tol,
        #                                          max.iter = max.iter,
        #                                          init = list(beta1 = b1.init.arr[,l,k],
        #                                                      beta2 = b2.init.arr[,l,k]),
        #                                          report.prog=TRUE)

        grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1.diag[-fold1],
                                                           rX1 = X1[-fold1,],
                                                           groups1 = groups1,
                                                           rY2 = Y2.diag[-fold2],
                                                           rX2 = X2[-fold2,],
                                                           groups2 = groups2,
                                                           rho1 = rho1,
                                                           rho2 = rho2,
                                                           lambda = lambda.seq[l]*(n.folds - 1)/n.folds,
                                                           eta = eta.seq[k]*(n.folds - 1)/n.folds,
                                                           w1 = w1,
                                                           w2 = w2,
                                                           w = w,
                                                           rAA1 = AA1,
                                                           rAA2 = AA2,
                                                           rCom = Com,
                                                           tol = tol,
                                                           maxiter = max.iter,
                                                           beta1_init = b1.init.arr[,l,k],
                                                           beta2_init = b2.init.arr[,l,k])

        b1.fold <- grouplasso2pop_logreg.out$beta1.hat
        b2.fold <- grouplasso2pop_logreg.out$beta2.hat

        b1.folds.arr[,l,k,fold] <- b1.fold
        b2.folds.arr[,l,k,fold] <- b2.fold

        iterations[step,2+fold] <- grouplasso2pop_logreg.out$iter

        P1.fold <- logit(X1[fold1,] %*% b1.fold)
        P2.fold <- logit(X2[fold2,] %*% b2.fold)

        minus2ll1.fold <- - 2 * rho1 * mean( Y1.diag[fold1] * log( P1.fold) + (1-Y1.diag[fold1]) * log( 1 - P1.fold))
        minus2ll2.fold <- - 2 * rho2 * mean( Y2.diag[fold2] * log( P2.fold) + (1-Y2.diag[fold2]) * log( 1 - P2.fold))

        minus2ll.arr[l,k,fold] <- minus2ll1.fold + minus2ll2.fold

      }

      print(c(l,k))
      step <- step + 1

    }

  }

  meanCVll <- apply(minus2ll.arr,c(1,2),mean)
  minimizers <- which(meanCVll == min(meanCVll), arr.ind=TRUE)

  which.lambda.cv <- minimizers[1]
  which.eta.cv <- minimizers[2]
  which.lambda.cv.under.zero.eta <- which.min(meanCVll[,1])

  output <- list( b1.folds.arr = b1.folds.arr,
                  b2.folds.arr = b2.folds.arr,
                  minus2ll.arr = minus2ll.arr,
                  which.lambda.cv = which.lambda.cv,
                  which.eta.cv = which.eta.cv,
                  lambda.seq = lambda.seq,
                  which.lambda.cv.under.zero.eta = which.lambda.cv.under.zero.eta,
                  eta.seq = eta.seq,
                  iterations = iterations)

  return(output)
}




#' Choose tuning parameters for two population group lasso estimator with group testing data
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values
#' @param n.eta the number of eta values
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#' 
#' 
#' @examples
#' grouplassogt2pop_data <- get_grouplassogt2pop_data(n1=400,n2=600)
#' 
#' grouplassogt2pop_cv.out <- grouplassogt2pop_cv(Y1 = grouplassogt2pop_data$Y1,
#'                                                Z1 = grouplassogt2pop_data$Z1,
#'                                                Se1 = grouplassogt2pop_data$Se1,
#'                                                Sp1 = grouplassogt2pop_data$Sp1,
#'                                                X1 = grouplassogt2pop_data$X1,
#'                                                groups1 = grouplassogt2pop_data$groups1,
#'                                                Y2 = grouplassogt2pop_data$Y2,
#'                                                Z2 = grouplassogt2pop_data$Z2,
#'                                                Se2 = grouplassogt2pop_data$Se2,
#'                                                Sp2 = grouplassogt2pop_data$Sp2,
#'                                                X2 = grouplassogt2pop_data$X2,
#'                                                groups2 = grouplassogt2pop_data$groups2,
#'                                                rho1 = 2,
#'                                                rho2 = 1,
#'                                                n.lambda = 5,
#'                                                n.eta = 5,
#'                                                lambda.min.ratio = 0.01,
#'                                                n.folds = 5,
#'                                                w1 = grouplassogt2pop_data$w1,
#'                                                w2 = grouplassogt2pop_data$w2,
#'                                                w = grouplassogt2pop_data$w,
#'                                                AA1 = grouplassogt2pop_data$AA1,
#'                                                AA2 = grouplassogt2pop_data$AA2,
#'                                                Com = grouplassogt2pop_data$Com,
#'                                                tol = 1e-2,
#'                                                max.iter = 500,
#'                                                report.prog = TRUE)
#' @export
grouplassogt2pop_cv <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,n.folds,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE){

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set,
  # which will be used as initial values for the crossvalidation training fits.
  grouplassogt2pop_grid.out <- grouplassogt2pop_grid(Y1 = Y1,
                                                     Z1 = Z1,
                                                     Se1 = Se1,
                                                     Sp1 = Sp1,
                                                     X1 = X1,
                                                     groups1 = groups1,
                                                     Y2 = Y2,
                                                     Z2 = Z2,
                                                     Se2 = Se2,
                                                     Sp2 = Sp2,
                                                     X2 = X2,
                                                     groups2 = groups2,
                                                     rho1 = rho1,
                                                     rho2 = rho2,
                                                     n.lambda = n.lambda,
                                                     n.eta = n.eta,
                                                     lambda.min.ratio = lambda.min.ratio,
                                                     w1 = w1,
                                                     w2 = w2,
                                                     w = w,
                                                     AA1 = AA1,
                                                     AA2 = AA2,
                                                     Com = Com,
                                                     E.approx = E.approx,
                                                     tol = tol,
                                                     max.iter = max.iter,
                                                     report.prog = report.prog)

  lambda.seq <- grouplassogt2pop_grid.out$lambda.seq
  eta.seq <- grouplassogt2pop_grid.out$eta.seq
  b1.arr <- grouplassogt2pop_grid.out$b1.arr
  b2.arr <- grouplassogt2pop_grid.out$b2.arr

  # do the crossvalidation
  grouplassogt2pop_cv_fixedgrid.out <- grouplassogt2pop_cv_fixedgrid(Y1 = Y1,
                                                                     Z1 = Z1,
                                                                     Se1 = Se1,
                                                                     Sp1 = Sp1,
                                                                     X1 = X1,
                                                                     groups1 = groups1,
                                                                     Y2 = Y2,
                                                                     Z2 = Z2,
                                                                     Se2 = Se2,
                                                                     Sp2 = Sp2,
                                                                     X2 = X2,
                                                                     groups2 = groups2,
                                                                     rho1 = rho1,
                                                                     rho2 = rho2,
                                                                     lambda.seq = grouplassogt2pop_grid.out$lambda.seq,
                                                                     eta.seq = grouplassogt2pop_grid.out$eta.seq,
                                                                     n.folds = n.folds,
                                                                     b1.init.arr = grouplassogt2pop_grid.out$b1.arr,
                                                                     b2.init.arr = grouplassogt2pop_grid.out$b2.arr,
                                                                     w1 = w1,
                                                                     w2 = w2,
                                                                     w = w,
                                                                     AA1 = AA1,
                                                                     AA2 = AA2,
                                                                     Com = Com,
                                                                     E.approx = E.approx,
                                                                     tol = tol,
                                                                     max.iter = max.iter)

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  b1.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplassogt2pop_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = grouplassogt2pop_cv_fixedgrid.out$iterations)

  return(output)

}


#' Choose tuning parameters for two population group lasso estimator with group testing data
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values
#' @param n.eta the number of eta values
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#'
#' @examples
#' grouplassogt2pop_data <- get_grouplassogt2pop_data(n1 = 1000, n2 = 1200)
#'
#' grouplassogt2pop_cv_adapt.out <- grouplassogt2pop_cv_adapt(Y1 = grouplassogt2pop_data$Y1,
#'                                                            Z1 = grouplassogt2pop_data$Z1,
#'                                                            Se1 = grouplassogt2pop_data$Se1,
#'                                                            Sp1 = grouplassogt2pop_data$Sp1,
#'                                                            X1 = grouplassogt2pop_data$X1,
#'                                                            groups1 = grouplassogt2pop_data$groups1,
#'                                                            Y2 = grouplassogt2pop_data$Y2,
#'                                                            Z2 = grouplassogt2pop_data$Z2,
#'                                                            Se2 = grouplassogt2pop_data$Se2,
#'                                                            Sp2 = grouplassogt2pop_data$Sp2,
#'                                                            X2 = grouplassogt2pop_data$X2,
#'                                                            groups2 = grouplassogt2pop_data$groups2,
#'                                                            rho1 = 2,
#'                                                            rho2 = 1,
#'                                                            n.lambda = 5,
#'                                                            n.eta = 5,
#'                                                            lambda.min.ratio = 0.01,
#'                                                            n.folds = 5,
#'                                                            w1 = grouplassogt2pop_data$w1,
#'                                                            w2 = grouplassogt2pop_data$w2,
#'                                                            w = grouplassogt2pop_data$w,
#'                                                            AA1 = grouplassogt2pop_data$AA1,
#'                                                            AA2 = grouplassogt2pop_data$AA2,
#'                                                            Com = grouplassogt2pop_data$Com,
#'                                                            tol = 1e-2,
#'                                                            max.iter = 500,
#'                                                            report.prog = TRUE)
#' @export
grouplassogt2pop_cv_adapt <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,n.folds,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE){


  # pull the individual diagnoses to be treated as true disease statuses to find the sequence of lambda values.
  Y1.diag <- pull.diagnoses(Z1,Y1)
  Y2.diag <- pull.diagnoses(Z2,Y2)

  # find lambda.max and lambda.min
  q1 <- length(unique(groups1))
  q2 <- length(unique(groups2))

  n1 <- nrow(X1)
  n2 <- nrow(X2)

  norms1 <- numeric(q1)
  norms2 <- numeric(q2)
  for(j in 2:max(q1,q2))
  {

    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)

    if(j <= q1){
      
      norms1[j] <- sqrt(sum((t(X1[,ind1]) %*% (Y1.diag - mean(Y1.diag)))^2)) / ( w1[j] * n1 / rho1)
      
    }
    if(j <= q2){
      
      norms2[j] <- sqrt(sum((t(X2[,ind2]) %*% (Y2.diag - mean(Y2.diag)))^2)) / (w2[j] * n2 / rho1)
      
    }

  }

  # smallest value of lambda which sets all the non-intercept entries of beta1 and beta2 equal to zero.
  lambda.max <- 2 * max(norms1,norms2)
  lambda.initial.fit <- lambda.min.ratio * lambda.max

  # fit a grouplasso2pop with eta = 0 and lambda as lambda.min.ratio*lambda.max.
  grouplassogt2pop.out <- grouplassogt2pop(Y1 = Y1,
                                           Z1 = Z1,
                                           Se1 = Se1,
                                           Sp1 = Sp1,
                                           X1 = X1,
                                           groups1 = groups1,
                                           Y2 = Y2,
                                           Z2 = Z2,
                                           Se2 = Se2,
                                           Sp2 = Sp2,
                                           X2 = X2,
                                           groups2 = groups2,
                                           rho1 = rho1,
                                           rho2 = rho2,
                                           lambda = lambda.initial.fit,
                                           eta = 0,
                                           w1 = w1,
                                           w2 = w2,
                                           w = w,
                                           AA1 = AA1,
                                           AA2 = AA2,
                                           Com = Com,
                                           E.approx = E.approx,
                                           tol = tol,
                                           max.iter = max.iter,
                                           report.prog = FALSE)

  # now make new values of w1, w2, and w based on these.
  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    w1[j] <- min(w1[j]/sqrt( sum( grouplassogt2pop.out$beta1.hat[ind]^2 )),1e10) #replace Inf with 1e10
  }
  for(j in 1:q2)
  {
    ind <- which(groups2 == j)
    w2[j] <- min(w2[j]/sqrt( sum( grouplassogt2pop.out$beta2.hat[ind]^2 )),1e10)
  }
  for( j in Com)
  {
    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)
    w[j] <- min(w[j]/sum( (AA1[[j]] %*% grouplassogt2pop.out$beta1.hat[ind1] - AA2[[j]] %*% grouplassogt2pop.out$beta2.hat[ind2] )^2),1e10)
  }

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplassogt2pop_grid.out <- grouplassogt2pop_grid(Y1 = Y1,
                                                     Z1 = Z1,
                                                     Se1 = Se1,
                                                     Sp1 = Sp1,
                                                     X1 = X1,
                                                     groups1 = groups1,
                                                     Y2 = Y2,
                                                     Z2 = Z2,
                                                     Se2 = Se2,
                                                     Sp2 = Sp2,
                                                     X2 = X2,
                                                     groups2 = groups2,
                                                     rho1 = rho1,
                                                     rho2 = rho2,
                                                     n.lambda = n.lambda,
                                                     n.eta = n.eta,
                                                     lambda.min.ratio = lambda.min.ratio,
                                                     w1 = w1,
                                                     w2 = w2,
                                                     w = w,
                                                     AA1 = AA1,
                                                     AA2 = AA2,
                                                     Com = Com,
                                                     E.approx = E.approx,
                                                     tol = tol,
                                                     max.iter = max.iter,
                                                     report.prog = TRUE)



  # do the crossvalidation
  grouplassogt2pop_cv_fixedgrid.out <- grouplassogt2pop_cv_fixedgrid(Y1 = Y1,
                                                                     Z1 = Z1,
                                                                     Se1 = Se1,
                                                                     Sp1 = Sp1,
                                                                     X1 = X1,
                                                                     groups1 = groups1,
                                                                     Y2 = Y2,
                                                                     Z2 = Z2,
                                                                     Se2 = Se2,
                                                                     Sp2 = Sp2,
                                                                     X2 = X2,
                                                                     groups2 = groups2,
                                                                     rho1 = rho1,
                                                                     rho2 = rho2,
                                                                     lambda.seq = grouplassogt2pop_grid.out$lambda.seq,
                                                                     eta.seq = grouplassogt2pop_grid.out$eta.seq,
                                                                     n.folds = n.folds,
                                                                     b1.init.arr = grouplassogt2pop_grid.out$b1.arr,
                                                                     b2.init.arr = grouplassogt2pop_grid.out$b2.arr,
                                                                     w1 = w1,
                                                                     w2 = w2,
                                                                     w = w,
                                                                     AA1 = AA1,
                                                                     AA2 = AA2,
                                                                     Com = Com,
                                                                     E.approx = E.approx,
                                                                     tol = tol,
                                                                     max.iter = max.iter)

  output <- list( b1.arr = grouplassogt2pop_grid.out$b1.arr,
                  b2.arr = grouplassogt2pop_grid.out$b2.arr,
                  b1.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplassogt2pop_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = grouplassogt2pop_grid.out$lambda.seq,
                  eta.seq = grouplassogt2pop_grid.out$eta.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w1 = w1,
                  w2 = w2,
                  w = w,
                  iterations = grouplassogt2pop_cv_fixedgrid.out$iterations)

  return(output)

}

#' Choose tuning parameters for two population group lasso estimator with group testing data
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param groups1 a vector indicating to which group each covariate of data set 2 belongs
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values
#' @param n.eta the number of eta values
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#'
#' @examples
#'
#' @export
grouplassogt2pop_cv_adapt_fixedgrid <- function(Y1,Z1,Se1,Sp1,X1,groups1,Y2,Z2,Se2,Sp2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,lambda.initial.fit,n.folds,w1,w2,w,AA1,AA2,Com,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE){

  # fit a grouplasso2pop with eta = 0 and lambda as lambda.min.ratio*lambda.max.
  grouplassogt2pop.out <- grouplassogt2pop(Y1 = Y1,
                                           Z1 = Z1,
                                           Se1 = Se1,
                                           Sp1 = Sp1,
                                           X1 = X1,
                                           groups1 = groups1,
                                           Y2 = Y2,
                                           Z2 = Z2,
                                           Se2 = Se2,
                                           Sp2 = Sp2,
                                           X2 = X2,
                                           groups2 = groups2,
                                           rho1 = rho1,
                                           rho2 = rho2,
                                           lambda = lambda.initial.fit,
                                           eta = 0,
                                           w1 = w1,
                                           w2 = w2,
                                           w = w,
                                           AA1 = AA1,
                                           AA2 = AA2,
                                           Com = Com,
                                           E.approx = E.approx,
                                           tol = tol,
                                           max.iter = max.iter,
                                           report.prog = FALSE)

  # now make new values of w1, w2, and w based on these.

  q1 <- length(unique(groups1))
  q2 <- length(unique(groups2))

  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    w1[j] <- min(w1[j]/sqrt( sum( grouplassogt2pop.out$beta1.hat[ind]^2 )),1e10) #replace Inf with 1e10
  }
  for(j in 1:q2)
  {
    ind <- which(groups2 == j)
    w2[j] <- min(w2[j]/sqrt( sum( grouplassogt2pop.out$beta2.hat[ind]^2 )),1e10)
  }
  for( j in Com)
  {
    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)
    w[j] <- min(w[j]/sum( (AA1[[j]] %*% grouplassogt2pop.out$beta1.hat[ind1] - AA2[[j]] %*% grouplassogt2pop.out$beta2.hat[ind2] )^2),1e10)
  }

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplassogt2pop_fixedgrid.out <- grouplassogt2pop_fixedgrid(Y1 = Y1,
                                                               Z1 = Z1,
                                                               Se1 = Se1,
                                                               Sp1 = Sp1,
                                                               X1 = X1,
                                                               groups1 = groups1,
                                                               Y2 = Y2,
                                                               Z2 = Z2,
                                                               Se2 = Se2,
                                                               Sp2 = Sp2,
                                                               X2 = X2,
                                                               groups2 = groups2,
                                                               rho1 = rho1,
                                                               rho2 = rho2,
                                                               lambda.seq = lambda.seq,
                                                               eta.seq = eta.seq,
                                                               w1 = w1,
                                                               w2 = w2,
                                                               w = w,
                                                               AA1 = AA1,
                                                               AA2 = AA2,
                                                               Com = Com,
                                                               E.approx = E.approx,
                                                               tol = tol,
                                                               max.iter = max.iter,
                                                               report.prog = TRUE)

  # do the crossvalidation
  grouplassogt2pop_cv_fixedgrid.out <- grouplassogt2pop_cv_fixedgrid(Y1 = Y1,
                                                                     Z1 = Z1,
                                                                     Se1 = Se1,
                                                                     Sp1 = Sp1,
                                                                     X1 = X1,
                                                                     groups1 = groups1,
                                                                     Y2 = Y2,
                                                                     Z2 = Z2,
                                                                     Se2 = Se2,
                                                                     Sp2 = Sp2,
                                                                     X2 = X2,
                                                                     groups2 = groups2,
                                                                     rho1 = rho1,
                                                                     rho2 = rho2,
                                                                     lambda.seq = lambda.seq,
                                                                     eta.seq = eta.seq,
                                                                     n.folds = n.folds,
                                                                     b1.init.arr = grouplassogt2pop_fixedgrid.out$b1.arr,
                                                                     b2.init.arr = grouplassogt2pop_fixedgrid.out$b2.arr,
                                                                     w1 = w1,
                                                                     w2 = w2,
                                                                     w = w,
                                                                     AA1 = AA1,
                                                                     AA2 = AA2,
                                                                     Com = Com,
                                                                     E.approx = E.approx,
                                                                     tol = tol,
                                                                     max.iter = max.iter)

  # collect output
  output <- list( b1.arr = grouplassogt2pop_fixedgrid.out$b1.arr,
                  b2.arr = grouplassogt2pop_fixedgrid.out$b2.arr,
                  b1.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplassogt2pop_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplassogt2pop_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w1 = w1,
                  w2 = w2,
                  w = w,
                  iterations = grouplassogt2pop_cv_fixedgrid.out$iterations)

  return(output)

}


#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#' semipaddgt2pop_data <- get_semipaddgt2pop_data(n1 = 500, n2 = 604)
#'
#' semipaddgt2pop.out <- semipaddgt2pop(Y1 = semipaddgt2pop_data$Y1,
#'                                      Z1 = semipaddgt2pop_data$Z1,
#'                                      Se1 = semipaddgt2pop_data$Se1,
#'                                      Sp1 = semipaddgt2pop_data$Sp1,
#'                                      X1 = semipaddgt2pop_data$X1,
#'                                      nonparm1 = semipaddgt2pop_data$nonparm1,
#'                                      Y2 = semipaddgt2pop_data$Y2,
#'                                      Z2 = semipaddgt2pop_data$Z2,
#'                                      Se2 = semipaddgt2pop_data$Se2,
#'                                      Sp2 = semipaddgt2pop_data$Sp2,
#'                                      X2 = semipaddgt2pop_data$X2,
#'                                      nonparm2 = semipaddgt2pop_data$nonparm2,
#'                                      rho1 = 2,
#'                                      rho2 = 1,
#'                                      w1 = 1,
#'                                      w2 = 1,
#'                                      w = 1,
#'                                      nCom = 4,
#'                                      d1 = semipaddgt2pop_data$nonparm1 * 15,
#'                                      d2 = semipaddgt2pop_data$nonparm2 * 10,
#'                                      xi = 1,
#'                                      lambda.beta = 1,
#'                                      lambda.f = 1,
#'                                      eta.beta = 1,
#'                                      eta.f = 1,
#'                                      tol = 1e-2,
#'                                      max.iter = 500,
#'                                      report.prog = TRUE)
#'
#' plot_semipaddgt2pop(semipaddgt2pop.out,
#'                     true.functions=list(f1 = semipaddgt2pop_data$f1,
#'                                         f2 = semipaddgt2pop_data$f2,
#'                                         X1 = semipaddgt2pop_data$X1,
#'                                         X2 = semipaddgt2pop_data$X2)
#' )
#' @export
semipaddgt2pop <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,E.approx = FALSE,tol=1e-3,max.iter=500,report.prog=FALSE)
{

  # prepare input for grouplassogt2pop function
  grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          nCom = nCom,
                                                          d1 = d1,
                                                          d2 = d2,
                                                          xi = xi,
                                                          w1 = w1,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)

  # get group lasso estimates using the EM-algorithm
  grouplassogt2pop.out <- grouplassogt2pop(Y1 = Y1,
                                           Z1 = Z1,
                                           Se1 = Se1,
                                           Sp1 = Sp1,
                                           X1 = grouplasso2pop_inputs$DD1.tilde,
                                           groups1 = grouplasso2pop_inputs$groups1,
                                           Y2 = Y2,
                                           Z2 = Z2,
                                           Se2 = Se2,
                                           Sp2 = Sp2,
                                           X2 = grouplasso2pop_inputs$DD2.tilde,
                                           groups2 = grouplasso2pop_inputs$groups2,
                                           rho1 = rho1,
                                           rho2 = rho2,
                                           lambda = grouplasso2pop_inputs$lambda,
                                           eta = grouplasso2pop_inputs$eta,
                                           w1 = grouplasso2pop_inputs$w1,
                                           w2 = grouplasso2pop_inputs$w2,
                                           w = grouplasso2pop_inputs$w,
                                           AA1 = grouplasso2pop_inputs$AA1.tilde,
                                           AA2 = grouplasso2pop_inputs$AA2.tilde,
                                           Com = grouplasso2pop_inputs$Com,
                                           E.approx = E.approx,
                                           tol = tol,
                                           max.iter = max.iter,
                                           init = NULL,
                                           report.prog = report.prog)

  # construct fitted functions from grouplasso2pop output
  semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          groups1 = grouplasso2pop_inputs$groups1,
                                                          knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                          emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                          QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                          b1 = grouplassogt2pop.out$beta1.hat,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          groups2 = grouplasso2pop_inputs$groups2,
                                                          knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                          emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                          QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                          b2 = grouplassogt2pop.out$beta2.hat)


  # collect output
  output <- list( f1.hat = semipaddgt2pop_fitted$f1.hat,
                  f2.hat = semipaddgt2pop_fitted$f2.hat,
                  f1.hat.design = semipaddgt2pop_fitted$f1.hat.design,
                  f2.hat.design = semipaddgt2pop_fitted$f2.hat.design,
                  rho1 = rho1,
                  rho2 = rho2,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  d1 = d1,
                  d2 = d2,
                  xi = xi,
                  knots.list1 = grouplasso2pop_inputs$knots.list1,
                  knots.list2 = grouplasso2pop_inputs$knots.list2,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  eta.beta = eta.beta,
                  eta.f = eta.f,
                  Com = grouplasso2pop_inputs$Com,
                  inner.iter = grouplassogt2pop.out$inner.iter)

  return(output)

}



#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param n.eta the number of eta values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#' semipaddgt2pop_data <- get_semipaddgt2pop_data(n1 = 500, n2 = 604)
#'
#' semipaddgt2pop_grid.out <- semipaddgt2pop_grid(Y1 = semipaddgt2pop_data$Y1,
#'                                                Z1 = semipaddgt2pop_data$Z1,
#'                                                Se1 = semipaddgt2pop_data$Se1,
#'                                                Sp1 = semipaddgt2pop_data$Sp1,
#'                                                X1 = semipaddgt2pop_data$X1,
#'                                                nonparm1 = semipaddgt2pop_data$nonparm1,
#'                                                Y2 = semipaddgt2pop_data$Y2,
#'                                                Z2 = semipaddgt2pop_data$Z2,
#'                                                Se2 = semipaddgt2pop_data$Se2,
#'                                                Sp2 = semipaddgt2pop_data$Sp2,
#'                                                X2 = semipaddgt2pop_data$X2,
#'                                                nonparm2 = semipaddgt2pop_data$nonparm2,
#'                                                rho1 = 2,
#'                                                rho2 = 1,
#'                                                w1 = 1,
#'                                                w2 = 1,
#'                                                w = 1,
#'                                                nCom = 4,
#'                                                d1 = semipaddgt2pop_data$nonparm1 * 15,
#'                                                d2 = semipaddgt2pop_data$nonparm2 * 10,
#'                                                xi = 1,
#'                                                n.lambda = 3,
#'                                                n.eta = 3,
#'                                                lambda.min.ratio=.01,
#'                                                lambda.beta = 1,
#'                                                lambda.f = 1,
#'                                                eta.beta = 1,
#'                                                eta.f = 1,
#'                                                tol = 1e-2,
#'                                                max.iter = 500,
#'                                                report.prog = TRUE)
#'
#' plot_semipaddgt2pop_grid(semipaddgt2pop_grid.out,
#'                          true.functions = list(f1 = semipaddgt2pop_data$f1,
#'                                                f2 = semipaddgt2pop_data$f2,
#'                                                X1 = semipaddgt2pop_data$X1,
#'                                                X2 = semipaddgt2pop_data$X2)
#' )
#' @export
semipaddgt2pop_grid <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplasso2pop function
  grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          nCom = nCom,
                                                          d1 = d1,
                                                          d2 = d2,
                                                          xi = xi,
                                                          w1 = w2,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)

  # get group lasso estimators over a grid of lambda and eta values
  grouplassogt2pop_grid.out <- grouplassogt2pop_grid(Y1 = Y1,
                                                     Z1 = Z1,
                                                     Se1 = Se1,
                                                     Sp1 = Sp1,
                                                     X1 = grouplasso2pop_inputs$DD1.tilde,
                                                     groups1 = grouplasso2pop_inputs$groups1,
                                                     Y2 = Y2,
                                                     Z2 = Z2,
                                                     Se2 = Se2,
                                                     Sp2 = Sp2,
                                                     X2 = grouplasso2pop_inputs$DD2.tilde,
                                                     groups2 =  grouplasso2pop_inputs$groups2,
                                                     rho1 = rho1,
                                                     rho2 = rho2,
                                                     n.lambda = n.lambda,
                                                     n.eta = n.eta,
                                                     lambda.min.ratio = lambda.min.ratio,
                                                     w1 = grouplasso2pop_inputs$w1,
                                                     w2 = grouplasso2pop_inputs$w2,
                                                     w = grouplasso2pop_inputs$w,
                                                     AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                     AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                     Com = grouplasso2pop_inputs$Com,
                                                     E.approx = E.approx,
                                                     tol = tol,
                                                     max.iter = max.iter,
                                                     report.prog = report.prog)

  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)
  for(l in 1:n.lambda)
  {
    f1.hat[[l]] <- vector("list",n.eta)
    f2.hat[[l]] <- vector("list",n.eta)
  }

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplassogt2pop_grid.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplassogt2pop_grid.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
                  f1.hat.design = f1.hat.design,
                  f2.hat.design = f2.hat.design,
                  rho1 = rho1,
                  rho2 = rho2,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  d1 = d1,
                  d2 = d2,
                  xi = xi,
                  knots.list1 = grouplasso2pop_inputs$knots.list1,
                  knots.list2 = grouplasso2pop_inputs$knots.list2,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  eta.beta = eta.beta,
                  eta.f = eta.f,
                  Com = grouplasso2pop_inputs$Com,
                  n.lambda = n.lambda,
                  n.eta = n.eta,
                  lambda.seq = grouplassogt2pop_grid.out$lambda.seq,
                  eta.seq = grouplassogt2pop_grid.out$eta.seq,
                  iterations = grouplassogt2pop_grid.out$iterations)

  return(output)

}

#' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param n.eta the number of eta values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#' semipaddgt2pop_data <- get_semipaddgt2pop_data(n1 = 500, n2 = 604)
#'
#' semipaddgt2pop_cv.out <- semipaddgt2pop_cv(Y1 = semipaddgt2pop_data$Y1,
#'                                            Z1 = semipaddgt2pop_data$Z1,
#'                                            Se1 = semipaddgt2pop_data$Se1,
#'                                            Sp1 = semipaddgt2pop_data$Sp1,
#'                                            X1 = semipaddgt2pop_data$X1,
#'                                            nonparm1 = semipaddgt2pop_data$nonparm1,
#'                                            Y2 = semipaddgt2pop_data$Y2,
#'                                            Z2 = semipaddgt2pop_data$Z2,
#'                                            Se2 = semipaddgt2pop_data$Se2,
#'                                            Sp2 = semipaddgt2pop_data$Sp2,
#'                                            X2 = semipaddgt2pop_data$X2,
#'                                            nonparm2 = semipaddgt2pop_data$nonparm2,
#'                                            rho1 = 2,
#'                                            rho2 = 1,
#'                                            w1 = 1,
#'                                            w2 = 1,
#'                                            w = 1,
#'                                            nCom = 4,
#'                                            d1 = semipaddgt2pop_data$nonparm1 * 15,
#'                                            d2 = semipaddgt2pop_data$nonparm2 * 10,
#'                                            xi = 1,
#'                                            n.lambda = 3,
#'                                            n.eta = 3,
#'                                            lambda.min.ratio =.001,
#'                                            lambda.beta = 1,
#'                                            lambda.f = 1,
#'                                            eta.beta = 1,
#'                                            eta.f = 1,
#'                                            tol = 1e-2,
#'                                            max.iter = 500,
#'                                            report.prog = TRUE)
#'
#' plot_semipaddgt2pop_cv(semipaddgt2pop_cv.out,
#'                        true.functions = list(f1 = semipaddgt2pop_data$f1,
#'                                              f2 = semipaddgt2pop_data$f2,
#'                                              X1 = semipaddgt2pop_data$X1,
#'                                              X2 = semipaddgt2pop_data$X2)
#' )
#' @export
semipaddgt2pop_cv <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplassogt2pop function
  grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          nCom = nCom,
                                                          d1 = d1,
                                                          d2 = d2,
                                                          xi = xi,
                                                          w1 = w1,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)

  # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
  grouplassogt2pop_cv.out <- grouplassogt2pop_cv(Y1 = Y1,
                                                 Z1 = Z1,
                                                 Se1 = Se1,
                                                 Sp1 = Sp1,
                                                 X1 = grouplasso2pop_inputs$DD1.tilde,
                                                 groups1 = grouplasso2pop_inputs$groups1,
                                                 Y2 = Y2,
                                                 Z2 = Z2,
                                                 Se2 = Se2,
                                                 Sp2 = Sp2,
                                                 X2 = grouplasso2pop_inputs$DD2.tilde,
                                                 groups2 = grouplasso2pop_inputs$groups2,
                                                 rho1 = rho1,
                                                 rho2 = rho2,
                                                 n.lambda = n.lambda,
                                                 n.eta = n.eta,
                                                 lambda.min.ratio = lambda.min.ratio,
                                                 n.folds = n.folds,
                                                 w1 = grouplasso2pop_inputs$w1,
                                                 w2 = grouplasso2pop_inputs$w2,
                                                 w = grouplasso2pop_inputs$w,
                                                 AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                 AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                 Com = grouplasso2pop_inputs$Com,
                                                 E.approx = E.approx,
                                                 tol = tol,
                                                 max.iter = max.iter,
                                                 report.prog = report.prog)

  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  f1.hat.folds <- vector("list",n.lambda)
  f2.hat.folds <- vector("list",n.lambda)

  for(l in 1:n.lambda)
  {
    f1.hat[[l]] <- vector("list",n.eta)
    f2.hat[[l]] <- vector("list",n.eta)

    f1.hat.folds[[l]] <- vector("list",n.eta)
    f2.hat.folds[[l]] <- vector("list",n.eta)

    for( k in 1:n.eta)
    {

      f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
      f2.hat.folds[[l]][[k]] <- vector("list",n.folds)

    }

  }

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplassogt2pop_cv.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplassogt2pop_cv.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

      for( fold in 1:n.folds)
      {

        semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                                nonparm1 = nonparm1,
                                                                groups1 = grouplasso2pop_inputs$groups1,
                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                                b1 = grouplassogt2pop_cv.out$b1.folds.arr[,l,k,fold],
                                                                X2 = X2,
                                                                nonparm2 = nonparm2,
                                                                groups2 = grouplasso2pop_inputs$groups2,
                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                                b2 = grouplassogt2pop_cv.out$b2.folds.arr[,l,k,fold])

        f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
        f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat

      }

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
                  f1.hat.folds = f1.hat.folds,
                  f2.hat.folds = f2.hat.folds,
                  f1.hat.design = f1.hat.design,
                  f2.hat.design = f2.hat.design,
                  beta1.hat = beta1.hat,
                  beta2.hat = beta2.hat,
                  rho1 = rho1,
                  rho2 = rho2,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  d1 = d1,
                  d2 = d2,
                  xi = xi,
                  knots.list1 = grouplasso2pop_inputs$knots.list1,
                  knots.list2 = grouplasso2pop_inputs$knots.list2,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  eta.beta = eta.beta,
                  eta.f = eta.f,
                  Com = grouplasso2pop_inputs$Com,
                  n.lambda = n.lambda,
                  n.eta = n.eta,
                  n.folds = n.folds,
                  lambda.seq = grouplassogt2pop_cv.out$lambda.seq,
                  eta.seq = grouplassogt2pop_cv.out$eta.seq,
                  which.lambda.cv = grouplassogt2pop_cv.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv.out$which.lambda.cv.under.zero.eta,
                  iterations = grouplassogt2pop_cv.out$iterations)

  class(output) <- "semipaddgt2pop_cv"

  return(output)

}



#' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param n.eta the number of eta values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#' semipaddgt2pop_data <- get_semipaddgt2pop_data(n1 = 500, n2 = 604)
#'
#' semipaddgt2pop_cv.out <- semipaddgt2pop_cv(Y1 = semipaddgt2pop_data$Y1,
#'                                            Z1 = semipaddgt2pop_data$Z1,
#'                                            Se1 = semipaddgt2pop_data$Se1,
#'                                            Sp1 = semipaddgt2pop_data$Sp1,
#'                                            X1 = semipaddgt2pop_data$X1,
#'                                            nonparm1 = semipaddgt2pop_data$nonparm1,
#'                                            Y2 = semipaddgt2pop_data$Y2,
#'                                            Z2 = semipaddgt2pop_data$Z2,
#'                                            Se2 = semipaddgt2pop_data$Se2,
#'                                            Sp2 = semipaddgt2pop_data$Sp2,
#'                                            X2 = semipaddgt2pop_data$X2,
#'                                            nonparm2 = semipaddgt2pop_data$nonparm2,
#'                                            rho1 = 2,
#'                                            rho2 = 1,
#'                                            w1 = 1,
#'                                            w2 = 1,
#'                                            w = 1,
#'                                            nCom = 4,
#'                                            d1 = semipaddgt2pop_data$nonparm1 * 15,
#'                                            d2 = semipaddgt2pop_data$nonparm2 * 10,
#'                                            xi = 1,
#'                                            n.lambda = 3,
#'                                            n.eta = 3,
#'                                            lambda.min.ratio =.001,
#'                                            lambda.beta = 1,
#'                                            lambda.f = 1,
#'                                            eta.beta = 1,
#'                                            eta.f = 1,
#'                                            tol = 1e-2,
#'                                            max.iter = 500,
#'                                            report.prog = TRUE)
#'
#' plot_semipaddgt2pop_cv(semipaddgt2pop_cv.out,
#'                        true.functions = list(f1 = semipaddgt2pop_data$f1,
#'                                              f2 = semipaddgt2pop_data$f2,
#'                                              X1 = semipaddgt2pop_data$X1,
#'                                              X2 = semipaddgt2pop_data$X2)
#' )
#' @export
semipaddgt2pop_cv_adapt <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplassogt2pop function
  grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          nCom = nCom,
                                                          d1 = d1,
                                                          d2 = d2,
                                                          xi = xi,
                                                          w1 = w1,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)

  # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
  grouplassogt2pop_cv_adapt.out <- grouplassogt2pop_cv_adapt(Y1 = Y1,
                                                             Z1 = Z1,
                                                             Se1 = Se1,
                                                             Sp1 = Sp1,
                                                             X1 = grouplasso2pop_inputs$DD1.tilde,
                                                             groups1 = grouplasso2pop_inputs$groups1,
                                                             Y2 = Y2,
                                                             Z2 = Z2,
                                                             Se2 = Se2,
                                                             Sp2 = Sp2,
                                                             X2 = grouplasso2pop_inputs$DD2.tilde,
                                                             groups2 = grouplasso2pop_inputs$groups2,
                                                             rho1 = rho1,
                                                             rho2 = rho2,
                                                             n.lambda = n.lambda,
                                                             n.eta = n.eta,
                                                             lambda.min.ratio = lambda.min.ratio,
                                                             n.folds = n.folds,
                                                             w1 = grouplasso2pop_inputs$w1,
                                                             w2 = grouplasso2pop_inputs$w2,
                                                             w = grouplasso2pop_inputs$w,
                                                             AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                             AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                             Com = grouplasso2pop_inputs$Com,
                                                             E.approx = E.approx,
                                                             tol = tol,
                                                             max.iter = max.iter,
                                                             report.prog = report.prog)

  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  f1.hat.folds <- vector("list",n.lambda)
  f2.hat.folds <- vector("list",n.lambda)

  for(l in 1:n.lambda)
  {
    f1.hat[[l]] <- vector("list",n.eta)
    f2.hat[[l]] <- vector("list",n.eta)

    f1.hat.folds[[l]] <- vector("list",n.eta)
    f2.hat.folds[[l]] <- vector("list",n.eta)

    for( k in 1:n.eta)
    {

      f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
      f2.hat.folds[[l]][[k]] <- vector("list",n.folds)

    }

  }

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplassogt2pop_cv_adapt.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplassogt2pop_cv_adapt.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

      for( fold in 1:n.folds)
      {

        semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                                nonparm1 = nonparm1,
                                                                groups1 = grouplasso2pop_inputs$groups1,
                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                                b1 = grouplassogt2pop_cv_adapt.out$b1.folds.arr[,l,k,fold],
                                                                X2 = X2,
                                                                nonparm2 = nonparm2,
                                                                groups2 = grouplasso2pop_inputs$groups2,
                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                                b2 = grouplassogt2pop_cv_adapt.out$b2.folds.arr[,l,k,fold])

        f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
        f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat

      }

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
                  f1.hat.folds = f1.hat.folds,
                  f2.hat.folds = f2.hat.folds,
                  f1.hat.design = f1.hat.design,
                  f2.hat.design = f2.hat.design,
                  beta1.hat = beta1.hat,
                  beta2.hat = beta2.hat,
                  rho1 = rho1,
                  rho2 = rho2,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  d1 = d1,
                  d2 = d2,
                  xi = xi,
                  knots.list1 = grouplasso2pop_inputs$knots.list1,
                  knots.list2 = grouplasso2pop_inputs$knots.list2,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  eta.beta = eta.beta,
                  eta.f = eta.f,
                  Com = grouplasso2pop_inputs$Com,
                  n.lambda = n.lambda,
                  n.eta = n.eta,
                  n.folds = n.folds,
                  lambda.seq = grouplassogt2pop_cv_adapt.out$lambda.seq,
                  eta.seq = grouplassogt2pop_cv_adapt.out$eta.seq,
                  lambda.initial.fit = grouplassogt2pop_cv_adapt.out$lambda.initial.fit,
                  which.lambda.cv = grouplassogt2pop_cv_adapt.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv_adapt.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv_adapt.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplassogt2pop_cv_adapt.out$w1,
                  w2 = grouplassogt2pop_cv_adapt.out$w2,
                  w = grouplassogt2pop_cv_adapt.out$w,
                  iterations = grouplassogt2pop_cv_adapt.out$iterations)

  class(output) <- "semipaddgt2pop_cv"

  return(output)

}



#' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se2 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp2 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param n.eta the number of eta values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#'
#' @export
semipaddgt2pop_cv_adapt_fixedgrid <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.seq,eta.seq,lambda.initial.fit,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplassogt2pop function
  grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
                                                          nonparm1 = nonparm1,
                                                          X2 = X2,
                                                          nonparm2 = nonparm2,
                                                          nCom = nCom,
                                                          d1 = d1,
                                                          d2 = d2,
                                                          xi = xi,
                                                          w1 = w1,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)

  # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
  grouplassogt2pop_cv_adapt_fixedgrid.out <- grouplassogt2pop_cv_adapt_fixedgrid(Y1 = Y1,
                                                                                 Z1 = Z1,
                                                                                 Se1 = Se1,
                                                                                 Sp1 = Sp1,
                                                                                 X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                                 groups1 = grouplasso2pop_inputs$groups1,
                                                                                 Y2 = Y2,
                                                                                 Z2 = Z2,
                                                                                 Se2 = Se2,
                                                                                 Sp2 = Sp2,
                                                                                 X2 = grouplasso2pop_inputs$DD2.tilde,
                                                                                 groups2 = grouplasso2pop_inputs$groups2,
                                                                                 rho1 = rho1,
                                                                                 rho2 = rho2,
                                                                                 lambda.seq = lambda.seq,
                                                                                 eta.seq = eta.seq,
                                                                                 lambda.initial.fit = lambda.initial.fit,
                                                                                 n.folds = n.folds,
                                                                                 w1 = grouplasso2pop_inputs$w1,
                                                                                 w2 = grouplasso2pop_inputs$w2,
                                                                                 w = grouplasso2pop_inputs$w,
                                                                                 AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                                                 AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                                                 Com = grouplasso2pop_inputs$Com,
                                                                                 E.approx = E.approx,
                                                                                 tol = tol,
                                                                                 max.iter = max.iter,
                                                                                 report.prog = report.prog)

  # get matrices of the fitted functions evaluated at the design points

  n.lambda <- length(lambda.seq)
  n.eta <- length(eta.seq)

  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  f1.hat.folds <- vector("list",n.lambda)
  f2.hat.folds <- vector("list",n.lambda)

  for(l in 1:n.lambda)
  {
    f1.hat[[l]] <- vector("list",n.eta)
    f2.hat[[l]] <- vector("list",n.eta)

    f1.hat.folds[[l]] <- vector("list",n.eta)
    f2.hat.folds[[l]] <- vector("list",n.eta)

    for( k in 1:n.eta)
    {

      f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
      f2.hat.folds[[l]][[k]] <- vector("list",n.folds)

    }

  }

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplassogt2pop_cv_adapt_fixedgrid.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplassogt2pop_cv_adapt_fixedgrid.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

      for( fold in 1:n.folds)
      {

        semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                                nonparm1 = nonparm1,
                                                                groups1 = grouplasso2pop_inputs$groups1,
                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                                b1 = grouplassogt2pop_cv_adapt_fixedgrid.out$b1.folds.arr[,l,k,fold],
                                                                X2 = X2,
                                                                nonparm2 = nonparm2,
                                                                groups2 = grouplasso2pop_inputs$groups2,
                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                                b2 = grouplassogt2pop_cv_adapt_fixedgrid.out$b2.folds.arr[,l,k,fold])

        f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
        f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat

      }

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
                  f1.hat.folds = f1.hat.folds,
                  f2.hat.folds = f2.hat.folds,
                  f1.hat.design = f1.hat.design,
                  f2.hat.design = f2.hat.design,
                  beta1.hat = beta1.hat,
                  beta2.hat = beta2.hat,
                  rho1 = rho1,
                  rho2 = rho2,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  d1 = d1,
                  d2 = d2,
                  xi = xi,
                  knots.list1 = grouplasso2pop_inputs$knots.list1,
                  knots.list2 = grouplasso2pop_inputs$knots.list2,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  eta.beta = eta.beta,
                  eta.f = eta.f,
                  Com = grouplasso2pop_inputs$Com,
                  n.lambda = n.lambda,
                  n.eta = n.eta,
                  n.folds = n.folds,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  which.lambda.cv = grouplassogt2pop_cv_adapt_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplassogt2pop_cv_adapt_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplassogt2pop_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplassogt2pop_cv_adapt_fixedgrid.out$w1,
                  w2 = grouplassogt2pop_cv_adapt_fixedgrid.out$w2,
                  w = grouplassogt2pop_cv_adapt_fixedgrid.out$w,
                  iterations = grouplassogt2pop_cv_adapt_fixedgrid.out$iterations)

  class(output) <- "semipaddgt2pop_cv"

  return(output)

}
