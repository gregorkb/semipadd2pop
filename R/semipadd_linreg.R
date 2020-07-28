#' Compute the objective function of the 1-population group lasso problem with a continuous response
#'
#' @param beta the vector of coefficients
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param lambda the level of sparsity penalization
#' @param w group-specific weights for different penalization of groups
#' @return Returns the value of the objective function for the group lasso problem
#' @export
grouplasso_linreg_obj <- function(beta,Y,X,groups,lambda,w)
{

  q <- length(unique(groups))
  n <- nrow(X)

  LScrit <- sum( (Y - X %*% beta)^2 )

  beta.wl2l1 <- 0
  for(j in 1:q)
  {
    ind <- which(groups == j)
    beta.wl2l1  <- beta.wl2l1 + w[j] * sqrt(sum( beta[ind]^2 ))
  }

  pen <- lambda * beta.wl2l1

  val <- LScrit + pen

  return(val)

}

#' Minimize the objective function of the 1-population group lasso problem with a continuous response
#'
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector of integers indicating to which group each covariate belongs
#' @param lambda the level of sparsity penalization
#' @param w group-specific weights for different penalization across groups
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the minimizer of the group lasso objective function
#'
#' @examples
#' # generate data
#' grouplasso_linreg_data <- get_grouplasso_linreg_data(n = 500)
#'
#' # fit grouplasso1pop_linreg estimator
#' grouplasso_linreg.out <- grouplasso_linreg_R(Y = grouplasso_linreg_data$Y,
#'                                              X = grouplasso_linreg_data$X,
#'                                              groups = grouplasso_linreg_data$groups,
#'                                              lambda = 10,
#'                                              w = grouplasso_linreg_data$w,
#'                                              tol = 1e-4,
#'                                              max.iter = 500,
#'                                              plot_obj = TRUE)
#' @export
grouplasso_linreg_R <- function(Y,X,groups,lambda,w,tol=1e-4,max.iter=500,plot_obj=FALSE,init = NULL)
{

  # initialize estimators and convergence criteria
  q <- ncol(X)
  n <- nrow(X)
  
  if( length(init) == 0){
    
    beta.hat1 <- rep(0,q)
    
  } else {
    
    beta.hat1 <- init
    
  }
  
  conv <- 1
  iter <- 0
  got.eigen <- numeric(q)
  
  eigen <- vector("list", q)
  
  
  d <- table(groups)
  singleton.grps <- which(d == 1)
  nonsingleton.grps <- which(d != 1)
  obj.val <- numeric()
  
  while( conv > tol & iter < max.iter)
  {
    
    beta.hat0 <- beta.hat1
    
    # go through the groups of size 1 of the first data set
    for(j in singleton.grps)
    {
      
      ind <- which(groups == j)
      rj <- as.numeric(Y - X[,-ind,drop=FALSE] %*% beta.hat1[-ind])
      sXj <- sum(X[,ind]^2)
      Xj.rj <- t(X[,ind]) %*% rj
      
      beta.hat1[ind] <- SoftThresh_R(Xj.rj,lambda*w[j])/sXj
        
    }
    
    # go through the groups of size more than 1 of the first data set
    for(j in nonsingleton.grps)
    {
      
      ind <- which(groups == j)
      rj <- as.numeric(Y - X[,-ind,drop=FALSE] %*% beta.hat1[-ind])
      Xrj.norm <- sqrt(sum(  (t(X[,ind]) %*% rj)^2 ))
        
        if( Xrj.norm < lambda * w[j])
        {
          
          beta.hat1[ind] <- 0
          
        } else {
          
          # compute FD update
          if(got.eigen[j]==0)
          {
            LtL <- t(X[,ind]) %*% X[,ind]
            eigen[[j]] <- eigen(LtL)
            got.eigen[j] <- 1
            
          }
          
          h <- rj
          L <- X[,ind]
          
          beta.hat1[ind] <- FoygelDrton_Armadillo(h,L,lambda*w[j],eigen[[j]]$values,t(eigen[[j]]$vectors))
          
        }
        
    } 
      
    conv <- max(abs(beta.hat1 - beta.hat0))
    iter <- iter + 1
    
    if(plot_obj == TRUE){
      
      obj.val[iter] <- grouplasso_linreg_obj(beta.hat1,Y,X,groups,lambda,w)
      
    }
    
  }
    
  if(plot_obj == TRUE){
    
    plot(obj.val)
    
  }
  
  beta.hat <- beta.hat1
  
  output <- list(beta.hat = beta.hat,
                 obj.val = obj.val,
                 iter = iter)
  
  return(output)

}

#' Fit group lasso regression estimator over a grid of lambda values
#'
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param w group-specific weights for different penalization of different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#' @examples
#' grouplasso_linreg_data <- get_grouplasso_linreg_data(n = 400)
#'
#' grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = grouplasso_linreg_data$Y,
#'                                                      X = grouplasso_linreg_data$X,
#'                                                      groups = grouplasso_linreg_data$groups,
#'                                                      n.lambda = 25,
#'                                                      lambda.min.ratio = 0.001,
#'                                                      lambda.max.ratio = 0.1,
#'                                                      w = grouplasso_linreg_data$w,
#'                                                      tol = 1e-3,
#'                                                      max.iter = 500,
#'                                                      report.prog = TRUE)
#' @export
grouplasso_linreg_grid <- function(Y,X,groups,n.lambda,lambda.min.ratio,lambda.max.ratio=1,w,tol=1e-4,max.iter=500,report.prog=FALSE)
{

  # find lambda.max
  q <- length(unique(groups))
  n <- nrow(X)

  norms <- numeric(q)
  for(j in 2:q)
  {

    ind <- which(groups == j)

    if(j <= q){

      norms[j] <- sqrt(sum((t(X[,ind]) %*% (Y - mean(Y)))^2)) / w[j]

    }

  }

  lambda.max <- 2 * max(norms) # yes this is correct! This is the smallest value of lambda which sets all the non-intercept entries of beta equal to zero.

  # make a lambda sequence
  largest.lambda <- lambda.max.ratio * lambda.max
  smallest.lambda <- lambda.min.ratio * lambda.max
  lambda.seq <- sort(c(exp(log(smallest.lambda) + ((n.lambda+1):1)/(n.lambda+1) * ((log(largest.lambda) - log(smallest.lambda)))))[-1])
  
  if(n.lambda == 1) lambda.seq <- lambda.min

  # fit over the values in the lambda sequence using warm starts
  b.mat <- matrix(0,ncol(X),n.lambda)
  iterations <- matrix(0,n.lambda,2)
  colnames(iterations) <- c("lambda","iter")
  step <- 0
  init <- rep(0,ncol(X))
  for(l in 1:n.lambda){

      grouplasso_linreg.out <- grouplasso_linreg(rY = Y,
                                                 rX = X,
                                                 groups = groups,
                                                 lambda = lambda.seq[l],
                                                 w = w,
                                                 tol = tol,
                                                 maxiter = max.iter,
                                                 beta_init = init)

      b <- grouplasso_linreg.out$beta.hat

      init <- b
      b.mat[,l] <- b

      step <- step + 1
      iterations[step,] <- c(lambda.seq[l],grouplasso_linreg.out$iter)

      if(report.prog == TRUE){

        print(c(l,grouplasso_linreg.out$iter))

      }

  }

  output <- list( b.mat = b.mat,
                  lambda.seq = lambda.seq,
                  iterations = iterations)

  return(output)

}

#' Choose tuning parameters by crossvalidation for grouplasso linreg when given a fixed grid of lambda values
#'
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param lambda.seq sequence of lambda values
#' @param n.folds the number of crossvalidation folds
#' @param b.init.arr array of initial values for beta
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization of different groups
#' @param tol the convergence tolerance
#' @param max.iter the maximum number of iterations allowed for each fit
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#' @examples
#' grouplasso_linreg_data <- get_grouplasso_linreg_data(n = 100)
#' 
#' grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = grouplasso_linreg_data$Y,
#'                                                      X = grouplasso_linreg_data$X,
#'                                                      groups = grouplasso_linreg_data$groups,
#'                                                      n.lambda = 25,
#'                                                      lambda.min.ratio = 0.001,
#'                                                      lambda.max.ratio = 0.1,
#'                                                      w = grouplasso_linreg_data$w,
#'                                                      tol = 1e-3,
#'                                                      max.iter = 500,
#'                                                      report.prog = FALSE)
#' 
#' grouplasso_linreg_cv_fixedgrid.out <- grouplasso_linreg_cv_fixedgrid(Y = grouplasso_linreg_data$Y,
#'                                                                      X = grouplasso_linreg_data$X,
#'                                                                      groups = grouplasso_linreg_data$groups,
#'                                                                      lambda.seq = grouplasso_linreg_grid.out$lambda.seq,
#'                                                                      n.folds = 5,
#'                                                                      b.init.mat = grouplasso_linreg_grid.out$b.mat,
#'                                                                      w = grouplasso_linreg_data$w,
#'                                                                      tol = 1e-3,
#'                                                                      max.iter = 500)
#' @export
grouplasso_linreg_cv_fixedgrid <- function(Y,X,groups,lambda.seq,n.folds,b.init.mat,w,tol=1e-3,max.iter=500)
{

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

        grouplasso_linreg.out <- grouplasso_linreg(rY = Y[-fold.ind],
                                                   rX = X[-fold.ind,],
                                                   groups = groups,
                                                   lambda = lambda.seq[l]*(n.folds - 1)/n.folds,
                                                   w = w,
                                                   tol = tol,
                                                   maxiter = max.iter,
                                                   beta_init = b.init.mat[,l])

        b.fold <- grouplasso_linreg.out$beta.hat
        b.folds.arr[,l,fold] <- b.fold

        minus2ll.fold <- sum( (Y[fold.ind] - X[fold.ind,] %*% b.fold )^2 )
        minus2ll.mat[l,fold] <- minus2ll.fold

        iterations[step,1+fold] <- grouplasso_linreg.out$iter

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


#' Choose tuning parameters by crossvalidation for grouplasso linreg.
#'
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#'
#' @examples
#' grouplasso_linreg_data <- get_grouplasso_linreg_data(n = 100)
#' 
#' grouplasso_linreg_cv.out <- grouplasso_linreg_cv(Y = grouplasso_linreg_data$Y,
#'                                                  X = grouplasso_linreg_data$X,
#'                                                  groups = grouplasso_linreg_data$groups,
#'                                                  n.lambda = 25,
#'                                                  lambda.min.ratio = 0.001,
#'                                                  lambda.max.ratio = 0.1,
#'                                                  n.folds = 5,
#'                                                  w = grouplasso_linreg_data$w,
#'                                                  tol = 1e-3,
#'                                                  max.iter = 500,
#'                                                  report.prog = FALSE)
#' @export
grouplasso_linreg_cv <- function(Y,X,groups,n.lambda,lambda.min.ratio,lambda.max.ratio=1,n.folds,w,tol=1e-4,max.iter=500,report.prog = TRUE){

  # obtain lambda.seq from the grid function, as well as the fits on the entire data set,
  # which will be used as initial values for the crossvalidation training fits.
  grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = Y,
                                                       X = X,
                                                       groups = groups,
                                                       n.lambda = n.lambda,
                                                       lambda.min.ratio = lambda.min.ratio,
                                                       lambda.max.ratio = lambda.max.ratio,
                                                       w = w,
                                                       tol = tol,
                                                       max.iter = max.iter,
                                                       report.prog = report.prog)

  lambda.seq <- grouplasso_linreg_grid.out$lambda.seq
  b.mat <- grouplasso_linreg_grid.out$b.mat

  # do the crossvalidation
  grouplasso_linreg_cv_fixedgrid.out <- grouplasso_linreg_cv_fixedgrid(Y = Y,
                                                                       X = X,
                                                                       groups = groups,
                                                                       lambda.seq = lambda.seq,
                                                                       n.folds = n.folds,
                                                                       b.init.mat = b.mat,
                                                                       w = w,
                                                                       tol = tol,
                                                                       max.iter = 500)

  output <- list( b.mat = b.mat,
                  b.folds.arr = grouplasso_linreg_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_linreg_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_linreg_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = lambda.seq,
                  iterations = grouplasso_linreg_cv_fixedgrid.out$iterations)

  return(output)

}


#' Choose tuning parameters by crossvalidation for grouplasso linreg with adaptive weights
#'
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#'
#' @examples
#' grouplasso_linreg_data <- get_grouplasso_linreg_data(n = 100)
#' 
#' grouplasso_linreg_cv_adapt.out <- grouplasso_linreg_cv_adapt(Y = grouplasso_linreg_data$Y,
#'                                                              X = grouplasso_linreg_data$X,
#'                                                              groups = grouplasso_linreg_data$groups,
#'                                                              n.lambda = 25,
#'                                                              lambda.min.ratio = 0.001,
#'                                                              lambda.max.ratio = 0.1,
#'                                                              n.folds = 5,
#'                                                              w = grouplasso_linreg_data$w,
#'                                                              tol = 1e-3,
#'                                                              max.iter = 500,
#'                                                              report.prog = FALSE)
#' @export
grouplasso_linreg_cv_adapt <- function(Y,X,groups,n.lambda,lambda.min.ratio,lambda.max.ratio=1,n.folds,w,tol=1e-3,max.iter=500,report.prog = TRUE){


  # find lambda.max
  q <- length(unique(groups))
  n <- nrow(X)

  norms <- numeric(q)
  for(j in 2:q){

    ind <- which(groups == j)

    if(j <= q){

      norms[j] <- sqrt(sum((t(X[,ind]) %*% (Y - mean(Y)))^2)) / w[j]

    }

  }

  lambda.max <- 2 * max(norms) # yes this is correct! This is the smallest value of lambda which sets all the non-intercept entries of beta equal to zero.
  lambda.initial.fit <- lambda.min.ratio * lambda.max

  # fit a grouplasso with lambda as lambda.min.ratio * lambda.max.
  grouplasso_linreg.out <- grouplasso_linreg(rY = Y,
                                             rX = X,
                                             groups = groups,
                                             lambda = lambda.initial.fit,
                                             w = w,
                                             tol = tol,
                                             maxiter = max.iter)


  # now make a new value of w based on this initial fit
  for(j in 1:q){

    ind <- which(groups == j)
    w[j] <- min(w[j]/sqrt( sum( grouplasso_linreg.out$beta.hat[ind]^2 )),1e10) # replace Inf with 1e10

  }

  # obtain lambda.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = Y,
                                                       X = X,
                                                       groups = groups,
                                                       n.lambda = n.lambda,
                                                       lambda.min.ratio = lambda.min.ratio,
                                                       lambda.max.ratio = lambda.max.ratio,
                                                       w = w,
                                                       tol = tol,
                                                       max.iter = max.iter,
                                                       report.prog = report.prog)

  lambda.seq <- grouplasso_linreg_grid.out$lambda.seq
  b.mat <- grouplasso_linreg_grid.out$b.mat

  # do the crossvalidation
  grouplasso_linreg_cv_fixedgrid.out <- grouplasso_linreg_cv_fixedgrid(Y = Y,
                                                                       X = X,
                                                                       groups = groups,
                                                                       lambda.seq = lambda.seq,
                                                                       n.folds = n.folds,
                                                                       b.init.mat = b.mat,
                                                                       w = w,
                                                                       tol = tol,
                                                                       max.iter = max.iter)

  output <- list( b.mat = b.mat,
                  b.folds.arr = grouplasso_linreg_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_linreg_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_linreg_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = lambda.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w = w,
                  iterations = grouplasso_linreg_cv_fixedgrid.out$iterations)

  return(output)

}

#' Fit semiparametric regression model on continuous-response data
#'
#' @param Y the response vector
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
#' semipadd_linreg_data <- get_semipadd_linreg_data(n = 200)
#' 
#' semipadd_linreg.out <- semipadd_linreg(Y = semipadd_linreg_data$Y,
#'                                        X = semipadd_linreg_data$X,
#'                                        nonparm = semipadd_linreg_data$nonparm,
#'                                        w = 1,
#'                                        d = 20,
#'                                        xi = 2,
#'                                        lambda.beta = 1,
#'                                        lambda.f = 1,
#'                                        tol = 1e-3,
#'                                        max.iter = 500,
#'                                        plot_obj = FALSE)
#' 
#' plot_semipaddgt(semipadd_linreg.out,
#'                 true.functions = list(f = semipadd_linreg_data$f,
#'                                       X = semipadd_linreg_data$X)
#'                 
#' )
#' @export
semipadd_linreg <- function(Y,X,nonparm,w,d,xi,lambda.beta,lambda.f,tol=1e-4,max.iter=500,plot_obj=FALSE)
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
  grouplasso_linreg.out <- grouplasso_linreg(rY = Y,
                                             rX = grouplasso_inputs$DD.tilde,
                                             groups = grouplasso_inputs$groups,
                                             lambda = grouplasso_inputs$lambda,
                                             w = grouplasso_inputs$w,
                                             tol = tol,
                                             maxiter = max.iter)

  # construct fitted functions from grouplasso output
  semipadd_fitted <- grouplasso_to_semipadd(X = X,
                                            nonparm = nonparm,
                                            groups = grouplasso_inputs$groups,
                                            knots.list = grouplasso_inputs$knots.list,
                                            emp.cent = grouplasso_inputs$emp.cent,
                                            QQ.inv = grouplasso_inputs$QQ.inv,
                                            b = grouplasso_linreg.out$beta.hat)

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


  class(output) <- "semipaddgt"

  return(output)

}
 


#' Compute semiparametric continuous-response regression model while penalizing over a grid of tuning parameter values
#'
#' @param Y the response vector
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
#' semipadd_linreg_data <- get_semipadd_linreg_data(n = 501)
#' 
#' semipadd_linreg_grid.out <- semipadd_linreg_grid(Y = semipadd_linreg_data$Y,
#'                                                  X = semipadd_linreg_data$X,
#'                                                  nonparm = semipadd_linreg_data$nonparm,
#'                                                  d = 15,
#'                                                  xi = 1,
#'                                                  w = 1,
#'                                                  lambda.beta = 1,
#'                                                  lambda.f = 1,
#'                                                  n.lambda = 10,
#'                                                  lambda.min.ratio = 0.001,
#'                                                  lambda.max.ratio = 0.5,
#'                                                  tol = 1e-3,
#'                                                  max.iter = 500,
#'                                                  report.prog = TRUE)
#' 
#' plot_semipaddgt_grid(semipadd_linreg_grid.out,
#'                      true.functions = list(f = semipadd_linreg_data$f,
#'                                            X = semipadd_linreg_data$X)
#'                      
#' )
#' @export
semipadd_linreg_grid <- function(Y,X,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,lambda.beta=1,lambda.f=1,tol=1e-3,max.iter = 1000,report.prog = FALSE)
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
  grouplasso_linreg_grid.out <- grouplasso_linreg_grid(Y = Y,
                                                       X = grouplasso_inputs$DD.tilde,
                                                       groups = grouplasso_inputs$groups,
                                                       n.lambda = n.lambda,
                                                       lambda.min.ratio = lambda.min.ratio,
                                                       lambda.max.ratio = lambda.max.ratio,
                                                       w = grouplasso_inputs$w,
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
                                                b = grouplasso_linreg_grid.out$b[,l])

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
                  lambda.seq = grouplasso_linreg_grid.out$lambda.seq,
                  iterations = grouplasso_linreg_grid.out$iterations)

  class(output) <- "semipaddgt_grid"

  return(output)

}


#' Compute semiparametric binary-response regression model while penalizing dissimilarity using CV to select tuning parameters
#'
#' @param Y the binary response vector
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
#' semipadd_linreg_data <- get_semipadd_linreg_data(n = 100)
#' 
#' semipadd_linreg_cv.out <- semipadd_linreg_cv(Y = semipadd_linreg_data$Y,
#'                                              X = semipadd_linreg_data$X,
#'                                              nonparm = semipadd_linreg_data$nonparm,
#'                                              d = 20,
#'                                              xi = 1,
#'                                              w = 1,
#'                                              lambda.beta = 1,
#'                                              lambda.f = 1,
#'                                              n.lambda = 25,
#'                                              lambda.min.ratio = 0.001,
#'                                              lambda.max.ratio = 0.01,
#'                                              tol = 1e-3,
#'                                              max.iter = 500,
#'                                              report.prog = TRUE)
#' 
#' plot_semipaddgt_cv(semipadd_linreg_cv.out,
#'                    true.functions = list(f = semipadd_linreg_data$f,
#'                                          X = semipadd_linreg_data$X)
#'                    
#' )
#' @export
semipadd_linreg_cv <- function(Y,X,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,max.iter = 1000,report.prog = FALSE)
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
  grouplasso_linreg_cv.out <- grouplasso_linreg_cv(Y = Y,
                                                   X = grouplasso_inputs$DD.tilde,
                                                   groups = grouplasso_inputs$groups,
                                                   n.lambda = n.lambda,
                                                   lambda.min.ratio = lambda.min.ratio,
                                                   lambda.max.ratio = lambda.max.ratio,
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
                                                  b = grouplasso_linreg_cv.out$b.mat[,l])

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
                                                    b = grouplasso_linreg_cv.out$b.folds.arr[,l,fold])

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
                  lambda.seq = grouplasso_linreg_cv.out$lambda.seq,
                  which.lambda.cv = grouplasso_linreg_cv.out$which.lambda.cv,
                  iterations = grouplasso_linreg_cv.out$iterations)

  class(output) <- "sempaddgt_cv"

  return(output)

}


#' Compute semiparametric continuous-response regression model using CV to select tuning parameters after an adaptive step
#'
#' @param Y the response vector
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
#' semipadd_linreg_data <- get_semipadd_linreg_data(n = 501)
#' 
#' semipadd_linreg_cv_adapt.out <- semipadd_linreg_cv_adapt(Y = semipadd_linreg_data$Y,
#'                                                          X = semipadd_linreg_data$X,
#'                                                          nonparm = semipadd_linreg_data$nonparm,
#'                                                          d = 25,
#'                                                          xi = 1,
#'                                                          w = 1,
#'                                                          lambda.beta = 1,
#'                                                          lambda.f = 1,
#'                                                          n.lambda = 25,
#'                                                          lambda.min.ratio = 0.001,
#'                                                          lambda.max.ratio = 0.1,
#'                                                          tol = 1e-3,
#'                                                          max.iter = 500,
#'                                                          report.prog = TRUE)
#' 
#' plot_semipaddgt_cv(semipadd_linreg_cv_adapt.out,
#'                    true.functions = list(f = semipadd_linreg_data$f,
#'                                          X = semipadd_linreg_data$X)
#'                    
#' )
#' @export
semipadd_linreg_cv_adapt <- function(Y,X,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio = 1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,max.iter = 1000,report.prog = FALSE)
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
  grouplasso_linreg_cv_adapt.out <- grouplasso_linreg_cv_adapt(Y = Y,
                                                               X = grouplasso_inputs$DD.tilde,
                                                               groups = grouplasso_inputs$groups,
                                                               n.lambda = n.lambda,
                                                               lambda.min.ratio = lambda.min.ratio,
                                                               lambda.max.ratio = lambda.max.ratio,
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
                                                b = grouplasso_linreg_cv_adapt.out$b.mat[,l])

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
                                                  b = grouplasso_linreg_cv_adapt.out$b.folds.arr[,l,fold])

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
                  lambda.seq = grouplasso_linreg_cv_adapt.out$lambda.seq,
                  which.lambda.cv = grouplasso_linreg_cv_adapt.out$which.lambda.cv,
                  lambda.initial.fit = grouplasso_linreg_cv_adapt.out$lambda.initial.fit,
                  iterations = grouplasso_linreg_cv_adapt.out$iterations)

  class(output) <- "sempaddgt_cv"

  return(output)

}