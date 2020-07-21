
#' Compute group lasso for two populations with group testing data
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param groups a vector indicating to which group each covariate belongs
#' @param lambda the level of sparsity penalization
#' @param w group-specific weights for different penalization of different groups
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param init a list of initial values for the coefficient 
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#' @examples 
#' grouplasso_gt_data <- get_grouplasso_gt_data(n = 1000)
#' 
#' grouplasso_gt.out <- grouplasso_gt(Y = grouplasso_gt_data$Y,
#'                                    Z = grouplasso_gt_data$Z,
#'                                    Se = grouplasso_gt_data$Se,
#'                                    Sp = grouplasso_gt_data$Sp,
#'                                    X = grouplasso_gt_data$X,
#'                                    groups = grouplasso_gt_data$groups,
#'                                    lambda = 1,
#'                                    w  = grouplasso_gt_data$w,
#'                                    E.approx = FALSE,
#'                                    tol = 1e-3,
#'                                    max.iter = 500,
#'                                    report.prog = TRUE)
#' @export
grouplasso_gt <- function(Y,Z,Se,Sp,X,groups,lambda,w,E.approx = FALSE,tol=1e-3,max.iter=1000,init=NULL,report.prog=FALSE)
{
  
  # set function to get conditional expectations
  get.EY <- eval(parse(text = ifelse(E.approx, "EYapprox","EYexact")))
  
  # set initial values
  if(length(init) == 0){
    
    beta.hat1 <- rep(0,ncol(X))
    
  } else{
    
    beta.hat1 <- init
    
  }
  
  ###### Do the EM-algorithm with penalized updates
  conv <- 1
  iter <- 0
  inner.iter <- numeric()
  while( conv > tol & iter < max.iter)
  {
    
    beta.hat0 <- beta.hat1
    
    
    # E-step: compute the conditional expectations for the true disease statuses
    EY <- as.numeric(get.EY(Z,Y,X=X,b=beta.hat1,Se,Sp))
    
    # update initial values
    init <- beta.hat1
    
    # M-step: maximize the objective function with conditional expectations substituted
    grouplasso_logreg.out <- grouplasso_logreg(rY = EY,
                                               rX = X,
                                               groups = groups,
                                               lambda = lambda,
                                               w = w,
                                               tol = tol,
                                               maxiter = max.iter,
                                               beta_init = beta.hat1)
    
    beta.hat1 <- grouplasso_logreg.out$beta.hat
    
    conv <- max(abs(beta.hat1 - beta.hat0))
    iter <- iter + 1
    inner.iter[iter] <- grouplasso_logreg.out$iter
    
    if(report.prog) print(grouplasso_logreg.out$iter)
    
  }
  
  output <- list( beta.hat = beta.hat1,
                  inner.iter = inner.iter)
  
}


#' Fit grouplasso logistic regression estimator with group testing data over a grid of lambda values
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param w group-specific weights for different penalization of different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#' @examples
#' grouplasso_gt_data <- get_grouplasso_gt_data(n = 1000)
#' 
#' grouplasso_gt_grid.out <- grouplasso_gt_grid(Y = grouplasso_gt_data$Y,
#'                                              Z = grouplasso_gt_data$Z,
#'                                              Se = grouplasso_gt_data$Se,
#'                                              Sp = grouplasso_gt_data$Sp,
#'                                              X = grouplasso_gt_data$X,
#'                                              groups = grouplasso_gt_data$groups,
#'                                              n.lambda = 10,
#'                                              lambda.min.ratio = 0.01,
#'                                              w  = grouplasso_gt_data$w,
#'                                              E.approx = FALSE,
#'                                              tol = 1e-3,
#'                                              max.iter = 500,
#'                                              report.prog = TRUE)
#' @export
grouplasso_gt_grid <- function(Y,Z,Se,Sp,X,groups,n.lambda,lambda.min.ratio,w,E.approx = FALSE,tol=1e-4,max.iter=500,report.prog=FALSE)
{
  # get diagnoses; determine lambda sequence using these.
  Y.diag <- pull.diagnoses(Z,Y)
  
  # find lambda.max
  q <- length(unique(groups))
  n <- nrow(X)

  norms <- numeric(q)
  for(j in 2:q)
  {

    ind <- which(groups == j)

    if(j <= q){

      norms[j] <- sqrt(sum((t(X[,ind]) %*% (Y.diag - mean(Y.diag)))^2)) / w[j]

    }

  }

  lambda.max <- 2 * max(norms) # yes this is correct! This is the smallest value of lambda which sets all the non-intercept entries of beta equal to zero.

  # make a lambda sequence
  lambda.min <- lambda.min.ratio * lambda.max
  lambda.seq <- sort(c(exp(log(lambda.min) + ((n.lambda+1):1)/(n.lambda+1) * ((log(lambda.max) - log(lambda.min)))))[-1])

  if(n.lambda == 1) lambda.seq <- lambda.min

  # fit over the values in the lambda sequence using warm starts
  b.mat <- matrix(0,ncol(X),n.lambda)
  iterations <- matrix(0,n.lambda,2)
  colnames(iterations) <- c("lambda","iter")
  step <- 0
  init <- NULL
  for(l in 1:n.lambda){

    grouplasso_logreg.out <- grouplasso_gt(Y = Y,
                                           Z = Z,
                                           Se = Se,
                                           Sp =Sp,
                                           X = X,
                                           groups = groups,
                                           lambda = lambda.seq[l],
                                           w = w,
                                           E.approx = E.approx,
                                           tol = tol,
                                           max.iter = max.iter,
                                           init = init)
                                               
    b <- grouplasso_logreg.out$beta.hat

    init <- b
    b.mat[,l] <- b

    step <- step + 1
    iterations[step,] <- c(lambda.seq[l],grouplasso_logreg.out$iter)

    if(report.prog == TRUE){

      print(c(l,grouplasso_logreg.out$iter))

    }

  }

  output <- list( b.mat = b.mat,
                  lambda.seq = lambda.seq,
                  iterations = iterations)

  return(output)

}



#' Choose tuning parameters by crossvalidation for grouplasso logreg when given a fixed grid of lambda values
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param lambda.seq sequence of lambda values
#' @param n.folds the number of crossvalidation folds
#' @param b.init.mat matrix of which the columns contain initial values for beta for the values of the tuning parameter in \code{lambda.seq}
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization of different groups
#' @param tol the convergence tolerance
#' @param max.iter the maximum number of iterations allowed for each fit
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#' @examples
#' grouplasso_gt_data <- get_grouplasso_gt_data(n = 400)
#' 
#' grouplasso_gt_grid.out <- grouplasso_gt_grid(Y = grouplasso_gt_data$Y,
#'                                              Z = grouplasso_gt_data$Z,
#'                                              Se = grouplasso_gt_data$Se,
#'                                              Sp = grouplasso_gt_data$Sp,
#'                                              X = grouplasso_gt_data$X,
#'                                              groups = grouplasso_gt_data$groups,
#'                                              n.lambda = 25,
#'                                              lambda.min.ratio = 0.01,
#'                                              w = grouplasso_gt_data$w,
#'                                              tol = 1e-3,
#'                                              max.iter = 500,
#'                                              report.prog = FALSE)
#' 
#' 
#' lambda.seq <- grouplasso_gt_grid.out$lambda.seq
#' b.mat <- grouplasso_gt_grid.out$b.mat
#' 
#' grouplasso_gt_cv_fixedgrid.out <- grouplasso_gt_cv_fixedgrid(Y = grouplasso_gt_data$Y,
#'                                                              Z = grouplasso_gt_data$Z,
#'                                                              Se = grouplasso_gt_data$Se,
#'                                                              Sp = grouplasso_gt_data$Sp,
#'                                                              X = grouplasso_gt_data$X,
#'                                                              groups = grouplasso_gt_data$groups,
#'                                                              lambda.seq = lambda.seq,
#'                                                              n.folds = 5,
#'                                                              b.init.mat = b.mat,
#'                                                              w = grouplasso_gt_data$w,
#'                                                              tol = 1e-3,
#'                                                              max.iter = 500)
#' @export
grouplasso_gt_cv_fixedgrid <- function(Y,Z,Se,Sp,X,groups,lambda.seq,n.folds,b.init.mat,w,tol=1e-3,max.iter=500)
{

  # for now just use the diagnoses as the true responses and use the logreg fits on these
  Y.diag <- pull.diagnoses(Z,Y)
  
  # create list of sets of indices indicating which observations are in each fold
  n <- nrow(X)

  folds <- vector("list",n.folds)
  fold.size <- floor(n / n.folds)
  for(fold in 1:n.folds){

    folds[[fold]] <- ((fold-1)*fold.size + 1):(fold*fold.size)

  }

  if( floor(n / n.folds) != n / n.folds ){

    folds[[n.folds]] <- c(folds[[n.folds]],(fold*fold.size+1):n)

  }

  # get fits at all lambda values on all cv folds
  n.lambda <- length(lambda.seq)
  b.folds.arr <- array(0,dim=c(ncol(X),n.lambda,n.folds))

  minus2ll.mat <- matrix(0,n.lambda,n.folds)

  iterations <- matrix(0,n.lambda,1+n.folds)
  colnames(iterations) <- c("lambda",paste("fold",1:n.folds,"iter"))
  step <- 1

  for(l in 1:n.lambda){

    iterations[step,1] <- lambda.seq[l]

    for(fold in 1:n.folds){

      fold.ind <- folds[[fold]]

      grouplasso_logreg.out <- grouplasso_logreg(rY = Y.diag[-fold.ind],
                                                 rX = X[-fold.ind,],
                                                 groups = groups,
                                                 lambda = lambda.seq[l]*(n.folds - 1)/n.folds,
                                                 w = w,
                                                 tol = tol,
                                                 maxiter = max.iter,
                                                 beta_init = b.init.mat[,l])
                                             

      b.fold <- grouplasso_logreg.out$beta.hat
      b.folds.arr[,l,fold] <- b.fold

      P.fold <- logit(X[fold.ind,] %*% b.fold)
      minus2ll.fold <- - 2 * sum( Y.diag[fold.ind] * log(P.fold) + (1-Y.diag[fold.ind]) * log( 1 - P.fold) )
      minus2ll.mat[l,fold] <- minus2ll.fold

      iterations[step,1+fold] <- grouplasso_logreg.out$iter

    }

    step <- step + 1

  }

  meanCVll <- apply(minus2ll.mat,1,mean)
  which.lambda.cv <- which.min(meanCVll)

  output <- list( b.folds.arr = b.folds.arr,
                  minus2ll.mat = minus2ll.mat,
                  which.lambda.cv = which.lambda.cv,
                  lambda.seq = lambda.seq,
                  iterations = iterations)

  return(output)

}


#' Choose tuning parameters for group lasso estimator with group testing data
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#' 
#' 
#' @examples
#' grouplasso_gt_data <- get_grouplasso_gt_data(n = 1000)
#' 
#' grouplasso_gt_cv.out <- grouplasso_gt_cv(Y = grouplasso_gt_data$Y,
#'                                          Z = grouplasso_gt_data$Z,
#'                                          Se = grouplasso_gt_data$Se,
#'                                          Sp = grouplasso_gt_data$Sp,
#'                                          X = grouplasso_gt_data$X,
#'                                          groups = grouplasso_gt_data$groups,
#'                                          n.lambda = 10,
#'                                          lambda.min.ratio = 0.01,
#'                                          n.folds = 5,
#'                                          w  = grouplasso_gt_data$w,
#'                                          E.approx = FALSE,
#'                                          tol = 1e-3,
#'                                          max.iter = 500,
#'                                          report.prog = TRUE)
#' @export
grouplasso_gt_cv <- function(Y,Z,Se,Sp,X,groups,n.lambda,lambda.min.ratio,n.folds,w,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE){
  
  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set,
  # which will be used as initial values for the crossvalidation training fits.
  grouplasso_gt_grid.out <- grouplasso_gt_grid(Y = Y,
                                             Z = Z,
                                             Se = Se,
                                             Sp = Sp,
                                             X = X,
                                             groups = groups,
                                             n.lambda = n.lambda,
                                             lambda.min.ratio = lambda.min.ratio,
                                             w = w,
                                             E.approx = E.approx,
                                             tol = tol,
                                             max.iter = max.iter,
                                             report.prog = report.prog)
  
  lambda.seq <- grouplasso_gt_grid.out$lambda.seq
  b.mat <- grouplasso_gt_grid.out$b.mat
  
  # do the crossvalidation
  grouplasso_gt_cv_fixedgrid.out <- grouplasso_gt_cv_fixedgrid(Y = Y,
                                                             Z = Z,
                                                             Se = Se,
                                                             Sp = Sp,
                                                             X = X,
                                                             groups = groups,
                                                             lambda.seq = grouplasso_gt_grid.out$lambda.seq,
                                                             n.folds = n.folds,
                                                             b.init.mat = grouplasso_gt_grid.out$b.mat,
                                                             w = w,
                                                             tol = tol,
                                                             max.iter = max.iter)
  
  output <- list( b.mat = b.mat,
                  b.folds.arr = grouplasso_gt_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_gt_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_gt_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = lambda.seq,
                  iterations = grouplasso_gt_cv_fixedgrid.out$iterations)
  
  return(output)
  
}


#' Choose tuning parameters for the group lasso estimator with group testing data
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the parametric model with group testing data
#'
#' @examples
#' grouplasso_gt_data <- get_grouplasso_gt_data(n = 1000)
#' 
#' grouplasso_gt_cv.out <- grouplasso_gt_cv_adapt(Y = grouplasso_gt_data$Y,
#'                                                Z = grouplasso_gt_data$Z,
#'                                                Se = grouplasso_gt_data$Se,
#'                                                Sp = grouplasso_gt_data$Sp,
#'                                                X = grouplasso_gt_data$X,
#'                                                groups = grouplasso_gt_data$groups,
#'                                                n.lambda = 10,
#'                                                lambda.min.ratio = 0.01,
#'                                                n.folds = 5,
#'                                                w  = grouplasso_gt_data$w,
#'                                                E.approx = FALSE,
#'                                                tol = 1e-3,
#'                                                max.iter = 500,
#'                                                report.prog = TRUE)
#' @export
grouplasso_gt_cv_adapt <- function(Y,Z,Se,Sp,X,groups,n.lambda,n.eta,lambda.min.ratio,n.folds,w,E.approx = FALSE,tol=1e-3,max.iter=1000,report.prog=TRUE){
  
  
  # pull the individual diagnoses to be treated as true disease statuses to find the sequence of lambda values.
  Y.diag <- pull.diagnoses(Z,Y)

  # find lambda.max and lambda.min
  q <- length(unique(groups))
  n <- nrow(X)

  norms <- numeric(q)
  for(j in 2:q)
  {
    
    ind <- which(groups == j)
    norms[j] <- sqrt(sum((t(X[,ind]) %*% (Y.diag - mean(Y.diag)))^2)) / w[j] 
    
  }
  
  # smallest value of lambda which sets all the non-intercept entries of beta1 and beta2 equal to zero.
  lambda.max <- 2 * max(norms)
  lambda.initial.fit <- lambda.min.ratio * lambda.max
  
  # fit a grouplasso_gt with eta = 0 and lambda as lambda.min.ratio*lambda.max.
  grouplasso_gt.out <- grouplasso_gt(Y = Y,
                                     Z = Z,
                                     Se = Se,
                                     Sp = Sp,
                                     X = X,
                                     groups = groups,
                                     lambda = lambda.initial.fit,
                                     w = w,
                                     E.approx = E.approx,
                                     tol = tol,
                                     max.iter = max.iter,
                                     report.prog = FALSE)
                                           
  # now make new values of w based on these.
  for(j in 1:q){
    
    ind <- which(groups == j)
    w[j] <- min(w[j]/sqrt( sum( grouplasso_gt.out$beta.hat[ind]^2 )),1e10) #replace Inf with 1e10
    
  }
  
  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplasso_gt_grid.out <- grouplasso_gt_grid(Y = Y,
                                               Z = Z,
                                               Se = Se,
                                               Sp = Sp,
                                               X = X,
                                               groups = groups,
                                               n.lambda = n.lambda,
                                               lambda.min.ratio = lambda.min.ratio,
                                               w = w,
                                               E.approx = E.approx,
                                               tol = tol,
                                               max.iter = max.iter,
                                               report.prog = report.prog)
  
  # do the crossvalidation
  grouplasso_gt_cv_fixedgrid.out <- grouplasso_gt_cv_fixedgrid(Y = Y,
                                                               Z = Z,
                                                               Se = Se,
                                                               Sp = Sp,
                                                               X = X,
                                                               groups = groups,
                                                               lambda.seq = grouplasso_gt_grid.out$lambda.seq,
                                                               n.folds = n.folds,
                                                               b.init.mat = grouplasso_gt_grid.out$b.mat,
                                                               w = w,
                                                               tol = tol,
                                                               max.iter = max.iter)
                                                                     
  output <- list( b.mat = grouplasso_gt_grid.out$b.mat,
                  b.folds.arr = grouplasso_gt_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_gt_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_gt_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = grouplasso_gt_grid.out$lambda.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w = w,
                  iterations = grouplasso_gt_cv_fixedgrid.out$iterations)
  
  return(output)
  
}


#' Compute semiparametric binary-response regression model with group-testing responses
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param d vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd_gt_data <- get_semipadd_gt_data(n = 500)
#' 
#' semipadd_gt.out <- semipadd_gt(Y = semipadd_gt_data$Y,
#'                                Z = semipadd_gt_data$Z,
#'                                Se = semipadd_gt_data$Se,
#'                                Sp = semipadd_gt_data$Sp,
#'                                X = semipadd_gt_data$X,
#'                                nonparm = semipadd_gt_data$nonparm,
#'                                w = 1,
#'                                d = semipadd_gt_data$nonparm*25,
#'                                xi = .5,
#'                                lambda.beta = 1,
#'                                lambda.f = 1,
#'                                tol = 1e-3,
#'                                max.iter = 500,
#'                                plot_obj = FALSE)
#' 
#' 
#' plot_semipaddgt(semipadd_gt.out,
#'                 true.functions = list(f = semipadd_gt_data$f,
#'                                       X = semipadd_gt_data$X)
#'                 
#' )
#' @export
semipadd_gt <- function(Y,Z,Se,Sp,X,nonparm,w,d,xi,lambda.beta,lambda.f,tol=1e-4,max.iter=500)
{

  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)

  # get group lasso estimators
  grouplasso_logreg.out <- grouplasso_gt(Y = Y,
                                         Z = Z,
                                         Se = Se,
                                         Sp = Sp,
                                         X = grouplasso_inputs$DD.tilde,
                                         groups = grouplasso_inputs$groups,
                                         lambda = grouplasso_inputs$lambda,
                                         w = grouplasso_inputs$w,
                                         tol = tol,
                                         max.iter = max.iter)
                                             

  # construct fitted functions from grouplasso output
  semipadd_fitted <- grouplasso_to_semipadd(X = X,
                                            nonparm = nonparm,
                                            groups = grouplasso_inputs$groups,
                                            knots.list = grouplasso_inputs$knots.list,
                                            emp.cent = grouplasso_inputs$emp.cent,
                                            QQ.inv = grouplasso_inputs$QQ.inv,
                                            b = grouplasso_logreg.out$beta.hat)

  # collect output
  output <- list(f.hat = semipadd_fitted$f.hat,
                 f.hat.design = semipadd_fitted$f.hat.design,
                 beta.hat = semipadd_fitted$beta.hat,
                 nonparm = nonparm,
                 d = d,
                 xi = xi,
                 knots.list = grouplasso_inputs$knots.list,
                 lambda.beta = lambda.beta,
                 lambda.f = lambda.f)


  class(output) <- "semipadd_gt"

  return(output)

}



#' Compute semiparametric binary-response regression model sets with group testing responses while penalizing over a grid of tuning parameter values
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization for different covariates
#' @param d vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd_gt_data <- get_semipadd_gt_data(n = 500)
#' 
#' semipadd_gt_grid.out <- semipadd_gt_grid(Y = semipadd_gt_data$Y,
#'                                          Z = semipadd_gt_data$Z,
#'                                          Se = semipadd_gt_data$Se,
#'                                          Sp = semipadd_gt_data$Sp,
#'                                          X = semipadd_gt_data$X,
#'                                          nonparm = semipadd_gt_data$nonparm,
#'                                          w = 1,
#'                                          d = semipadd_gt_data$nonparm*25,
#'                                          xi = 1,
#'                                          n.lambda = 20,
#'                                          lambda.min.ratio = .001,
#'                                          lambda.beta = 1,
#'                                          lambda.f = 1,
#'                                          tol = 1e-3,
#'                                          max.iter = 500,
#'                                          report.prog = FALSE)
#' 
#' plot_semipaddgt_grid(semipadd_gt_grid.out,
#'                      true.functions = list(f = semipadd_gt_data$f,
#'                                            X = semipadd_gt_data$X)
#'                      
#' )
#' @export
semipadd_gt_grid <- function(Y,Z,Se,Sp,X,nonparm,w,nCom,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.beta=1,lambda.f=1,E.approx = FALSE, tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)

  # get group lasso estimators over a grid of lambda and eta values
  grouplasso_gt_grid.out <- grouplasso_gt_grid(Y = Y,
                                               Z = Z,
                                               Se = Se,
                                               Sp = Se,
                                               X = grouplasso_inputs$DD.tilde,
                                               groups = grouplasso_inputs$groups,
                                               n.lambda = n.lambda,
                                               lambda.min.ratio = lambda.min.ratio,
                                               w = grouplasso_inputs$w,
                                               E.approx = E.approx,
                                               tol = tol,
                                               max.iter = max.iter,
                                               report.prog = report.prog)
                                                       
  # get matrices of the fitted functions evaluated at the design points
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X),n.lambda))
  beta.hat <- matrix(0,ncol(X),n.lambda)

  f.hat <- vector("list",n.lambda)

  for(l in 1:n.lambda){

    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_gt_grid.out$b[,l])

    f.hat[[l]] <- semipaddgt_fitted$f.hat
    f.hat.design[,,l] <- semipaddgt_fitted$f.hat.design
    beta.hat[,l] <- semipaddgt_fitted$beta.hat

  }

  # prepare output
  output <- list( f.hat = f.hat,
                  f.hat.design = f.hat.design,
                  beta.hat = beta.hat,
                  nonparm = nonparm,
                  d = d,
                  xi = xi,
                  knots.list = grouplasso_inputs$knots.list,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  n.lambda = n.lambda,
                  lambda.seq = grouplasso_gt_grid.out$lambda.seq,
                  iterations = grouplasso_gt_grid.out$iterations)

  class(output) <- "semipadd_gt_grid"

  return(output)

}



#' Compute semiparametric binary-response regression model with group testing responses using CV to select tuning parameters
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param d vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd_gt_data <- get_semipadd_gt_data(n = 500)
#' 
#' semipadd_gt_cv.out <- semipadd_gt_cv(Y = semipadd_gt_data$Y,
#'                                      Z = semipadd_gt_data$Z,
#'                                      Se = semipadd_gt_data$Se,
#'                                      Sp = semipadd_gt_data$Sp,
#'                                      X = semipadd_gt_data$X,
#'                                      nonparm = semipadd_gt_data$nonparm,
#'                                      w = 1,
#'                                      d = semipadd_gt_data$nonparm*25,
#'                                      xi = 1,
#'                                      n.lambda = 20,
#'                                      lambda.min.ratio = .001,
#'                                      n.folds = 5,
#'                                      lambda.beta = 1,
#'                                      lambda.f = 1,
#'                                      tol = 1e-3,
#'                                      max.iter = 500,
#'                                      report.prog = FALSE)
#' 
#' plot_semipaddgt_cv(semipadd_gt_cv.out,
#'                    true.functions = list(f = semipadd_gt_data$f,
#'                                          X = semipadd_gt_data$X)
#' )
#' @export
semipadd_gt_cv <- function(Y,Z,Se,Sp,X,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)

  # get group lasso estimators over a grid of lambda and eta values
  grouplasso_gt_cv.out <- grouplasso_gt_cv(Y = Y,
                                           Z = Z,
                                           Se = Se,
                                           Sp = Sp,
                                           X = grouplasso_inputs$DD.tilde,
                                           groups = grouplasso_inputs$groups,
                                           n.lambda = n.lambda,
                                           lambda.min.ratio = lambda.min.ratio,
                                           n.folds = n.folds,
                                           w = grouplasso_inputs$w,
                                           tol = tol,
                                           max.iter = max.iter,
                                           report.prog = report.prog)
                                                   
  # get matrices of the fitted functions evaluated at the design points
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X),n.lambda))
  beta.hat <- matrix(0,ncol(X),n.lambda)

  f.hat <- vector("list",n.lambda)
  f.hat.folds <- vector("list",n.lambda)

  for(l in 1:n.lambda){

    f.hat[[l]] <- vector("list",n.folds)
    f.hat.folds[[l]] <- vector("list",n.folds)

  }

  for(l in 1:n.lambda){

    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_gt_cv.out$b.mat[,l])

    f.hat[[l]] <- semipaddgt_fitted$f.hat
    f.hat.design[,,l] <- semipaddgt_fitted$f.hat.design
    beta.hat[,l] <- semipaddgt_fitted$beta.hat

    for( fold in 1:n.folds){

      semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                  nonparm = nonparm,
                                                  groups = grouplasso_inputs$groups,
                                                  knots.list = grouplasso_inputs$knots.list,
                                                  emp.cent = grouplasso_inputs$emp.cent,
                                                  QQ.inv = grouplasso_inputs$QQ.inv,
                                                  b = grouplasso_gt_cv.out$b.folds.arr[,l,fold])

      f.hat.folds[[l]][[fold]] <- semipaddgt_fitted$f.hat

    }

  }

  # prepare output
  output <- list( f.hat = f.hat,
                  f.hat.folds = f.hat.folds,
                  f.hat.design = f.hat.design,
                  beta.hat = beta.hat,
                  nonparm = nonparm,
                  d = d,
                  xi = xi,
                  knots.list = grouplasso_inputs$knots.list,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  n.lambda = n.lambda,
                  n.folds = n.folds,
                  lambda.seq = grouplasso_gt_cv.out$lambda.seq,
                  which.lambda.cv = grouplasso_gt_cv.out$which.lambda.cv,
                  iterations = grouplasso_gt_cv.out$iterations)

  class(output) <- "sempadd_gt_cv"

  return(output)

}


#' Compute semiparametric binary-response regression model with group testing responses using CV to select the tuning parameter after an adaptive step
#'
#' @param Y Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output in the format as output by one of the functions \code{individual.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param d vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd_gt_data <- get_semipadd_gt_data(n = 500)
#' 
#' semipadd_gt_cv.out <- semipadd_gt_cv(Y = semipadd_gt_data$Y,
#'                                      Z = semipadd_gt_data$Z,
#'                                      Se = semipadd_gt_data$Se,
#'                                      Sp = semipadd_gt_data$Sp,
#'                                      X = semipadd_gt_data$X,
#'                                      nonparm = semipadd_gt_data$nonparm,
#'                                      w = 1,
#'                                      d = semipadd_gt_data$nonparm*25,
#'                                      xi = 1,
#'                                      n.lambda = 20,
#'                                      lambda.min.ratio = .001,
#'                                      n.folds = 5,
#'                                      lambda.beta = 1,
#'                                      lambda.f = 1,
#'                                      tol = 1e-3,
#'                                      max.iter = 500,
#'                                      report.prog = FALSE)
#' 
#' 
#' plot_semipaddgt_cv(semipadd_gt_cv.out,
#'                    true.functions = list(f = semipadd_gt_data$f,
#'                                          X = semipadd_gt_data$X)
#'                    
#' )
#' @export
semipadd_gt_cv_adapt <- function(Y,Z,Se,Sp,X,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,max.iter = 1000,report.prog = FALSE)
{

  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)

  # get group lasso estimators over a grid of lambda and eta values
  grouplasso_logreg_cv_adapt.out <- grouplasso_gt_cv_adapt(Y = Y,
                                                           Z = Z,
                                                           Se = Se,
                                                           Sp = Sp,
                                                           X = grouplasso_inputs$DD.tilde,
                                                           groups = grouplasso_inputs$groups,
                                                           n.lambda = n.lambda,
                                                           lambda.min.ratio = lambda.min.ratio,
                                                           n.folds = n.folds,
                                                           w = grouplasso_inputs$w,
                                                           tol = tol,
                                                           max.iter = max.iter,
                                                           report.prog = report.prog)
                                                               
  # get matrices of the fitted functions evaluated at the design points
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X),n.lambda))
  beta.hat <- matrix(0,ncol(X),n.lambda)

  f.hat <- vector("list",n.lambda)
  f.hat.folds <- vector("list",n.lambda)

  for(l in 1:n.lambda){

    f.hat[[l]] <- vector("list",n.folds)
    f.hat.folds[[l]] <- vector("list",n.folds)

  }

  for(l in 1:n.lambda){

    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_logreg_cv_adapt.out$b.mat[,l])

    f.hat[[l]] <- semipaddgt_fitted$f.hat
    f.hat.design[,,l] <- semipaddgt_fitted$f.hat.design
    beta.hat[,l] <- semipaddgt_fitted$beta.hat

    for( fold in 1:n.folds){

      semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                  nonparm = nonparm,
                                                  groups = grouplasso_inputs$groups,
                                                  knots.list = grouplasso_inputs$knots.list,
                                                  emp.cent = grouplasso_inputs$emp.cent,
                                                  QQ.inv = grouplasso_inputs$QQ.inv,
                                                  b = grouplasso_logreg_cv_adapt.out$b.folds.arr[,l,fold])

      f.hat.folds[[l]][[fold]] <- semipaddgt_fitted$f.hat

    }

  }

  # prepare output
  output <- list( f.hat = f.hat,
                  f.hat.folds = f.hat.folds,
                  f.hat.design = f.hat.design,
                  beta.hat = beta.hat,
                  nonparm = nonparm,
                  d = d,
                  xi = xi,
                  knots.list = grouplasso_inputs$knots.list,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  n.lambda = n.lambda,
                  n.folds = n.folds,
                  lambda.seq = grouplasso_logreg_cv_adapt.out$lambda.seq,
                  which.lambda.cv = grouplasso_logreg_cv_adapt.out$which.lambda.cv,
                  lambda.initial.fit = grouplasso_logreg_cv_adapt.out$lambda.initial.fit,
                  iterations = grouplasso_logreg_cv_adapt.out$iterations)

  class(output) <- "sempadd_gt_cv"

  return(output)

}