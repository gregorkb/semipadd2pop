#' Compute the objective function of the 1-population group lasso problem with a binary response
#'
#' @param beta the vector of coefficients 
#' @param Y the response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate of data set 1 belongs
#' @param lambda the level of sparsity penalization
#' @param w group-specific weights for different penalization of groups
#' @return Returns the value of the objective function for the 2-population group lasso problem
#' @export
grouplasso_logreg_obj <- function(beta,Y,X,groups,lambda,w)
{

  q <- length(unique(groups))
  n <- nrow(X)

  P <- logit(X %*% beta)

  neg2LL1 <- - 2 * sum(  Y * log(P) + (1 - Y)*log(1-P))

  beta.wl2l1 <- 0
  for(j in 1:q)
  {
    ind <- which(groups == j)
    beta.wl2l1  <- beta.wl2l1 + w[j] * sqrt(sum( beta[ind]^2 ))
  }

  pen <- lambda * beta.wl2l1

  val <- neg2LL1 + pen

  return(val)

}
#' Minimize the objective function of the 1-population group lasso problem with a binary response
#'
#' @param Y the binary response vector
#' @param X matrix containing the design matrices
#' @param groups a vector of integers indicating to which group each covariate belongs
#' @param lambda the level of sparsity penalization
#' @param w1 group-specific weights for different penalization across groups
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the minimizer of the group lasso objective function
#'
#' @examples
#' # generate data
#' grouplasso_logreg_data <- get_grouplasso_logreg_data(n = 500)
#' 
#' # fit grouplasso1pop_logreg estimator
#' grouplasso_logreg.out <- grouplasso_logreg_R(Y = grouplasso_logreg_data$Y,
#'                                              X = grouplasso_logreg_data$X,
#'                                              groups = grouplasso_logreg_data$groups,
#'                                              lambda = 10,
#'                                              w = grouplasso_logreg_data$w,
#'                                              tol = 1e-4,
#'                                              max.iter = 500,
#'                                              plot_obj = TRUE)
#' @export
grouplasso_logreg_R <- function(Y,X,groups,lambda,w,tol=1e-4,max.iter=500,plot_obj=FALSE,init = NULL)
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

    # go through the groups of size 1
    for(j in singleton.grps)
    {
      
      P <- logit( X %*% beta.hat1)
      if( any(P  > 1 - 1e-10) | any(P < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # fixed Hessian
      ind <- which(groups == j)
      Xj <- X[,ind,drop=FALSE]
      Xj.tilde <- (1/2) * Xj
      Zj.tilde <- (1/2) * ( Xj %*% beta.hat1[ind] + 4 * ( Y - P ) )
      sxj <- sum(Xj.tilde^2)
      
      sxzj <- sum(Xj.tilde * Zj.tilde)
      beta.hat1[ind] <- SoftThresh_R(sxzj,lambda*w[j]/2)/sxj
      
    }

    # go through the groups of size more than 1
    for(j in nonsingleton.grps)
    {
      
      P <- logit( X %*% beta.hat1)
      if( any(P  > 1 - 1e-10) | any(P < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind <- which(groups == j)
      Xj <- X[,ind,drop=FALSE]
      
      # fixed Hessian approach for the groups
      Xj.tilde <- (1/2)*Xj
      Zj.tilde <- (1/2)*(Xj %*% beta.hat1[ind] + 4*(Y - P))
      xzj.norm <- sqrt(sum(  (t(Xj.tilde) %*% Zj.tilde)^2 ))
      
      if( xzj.norm < lambda * w[j]/2)
      {
        
        beta.hat1[ind] <- 0
        
      } else {
        
        # compute FD update
        if(got.eigen[j]==0){
          
          LtL <- t(Xj.tilde) %*% Xj.tilde
          eigen[[j]] <- eigen(LtL)
          got.eigen[j] <- 1
          
        }
        
        h <- Zj.tilde
        L <- Xj.tilde
        
        beta.hat1[ind] <- FoygelDrton_R(h,L,lambda*w[j]/2,eigen[[j]]$values,t(eigen[[j]]$vectors))
        
      }
      
    } 
    
    conv <- max(abs(c(beta.hat1 - beta.hat0)))
    iter <- iter + 1

    if(plot_obj == TRUE)
    {
      obj.val[iter] <- grouplasso_logreg_obj(beta.hat1,Y,X,groups,lambda,w)
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

#' Fit grouplasso logistic regression estimator over a grid of lambda and eta values
#'
#' @param Y the binary response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param w group-specific weights for different penalization of different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#' @examples 
#' grouplasso_logreg_data <- get_grouplasso_logreg_data(n = 400)
#' 
#' grouplasso_logreg_grid.out <- grouplasso_logreg_grid(Y = grouplasso_logreg_data$Y,
#'                                                      X = grouplasso_logreg_data$X,
#'                                                      groups = grouplasso_logreg_data$groups,
#'                                                      n.lambda = 25,
#'                                                      lambda.min.ratio = 0.001,
#'                                                      w = grouplasso_logreg_data$w,
#'                                                      tol = 1e-3,
#'                                                      max.iter = 500,
#'                                                      report.prog = TRUE)
#' @export
grouplasso_logreg_grid <- function(Y,X,groups,n.lambda,lambda.min.ratio,w,tol=1e-4,max.iter=500,report.prog=FALSE)
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
  lambda.min <- lambda.min.ratio * lambda.max
  lambda.seq <- sort(c(exp(log(lambda.min) + ((n.lambda+1):1)/(n.lambda+1) * ((log(lambda.max) - log(lambda.min)))))[-1])
  
  if(n.lambda == 1) lambda.seq <- lambda.min
  
  # fit over the values in the lambda sequence using warm starts
  b.mat <- matrix(0,ncol(X),n.lambda)
  iterations <- matrix(0,n.lambda,2)
  colnames(iterations) <- c("lambda","iter")
  step <- 0
  init <- rep(0,q)
  for(l in 1:n.lambda){
      
      grouplasso_logreg.out <- grouplasso_logreg(rY = Y,
                                                 rX = X,
                                                 groups = groups,
                                                 lambda = lambda.seq[l],
                                                 w = w,
                                                 tol = tol,
                                                 maxiter = max.iter,
                                                 beta_init = init)
                                                         
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
#' @param Y the binary response vector
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
#' grouplasso_logreg_data <- get_grouplasso_logreg_data(n = 100)
#' 
#' grouplasso_logreg_grid.out <- grouplasso_logreg_grid(Y = grouplasso_logreg_data$Y,
#'                                                      X = grouplasso_logreg_data$X,
#'                                                      groups = grouplasso_logreg_data$groups,
#'                                                      n.lambda = 25,
#'                                                      lambda.min.ratio = 0.01,
#'                                                      w = grouplasso_logreg_data$w,
#'                                                      tol = 1e-3,
#'                                                      max.iter = 500,
#'                                                      report.prog = FALSE)
#' 
#' lambda.seq <- grouplasso_logreg_grid.out$lambda.seq
#' b.mat <- grouplasso_logreg_grid.out$b.mat
#' 
#' grouplasso_logreg_cv_fixedgrid.out <- grouplasso_logreg_cv_fixedgrid(Y = grouplasso_logreg_data$Y,
#'                                                                      X = grouplasso_logreg_data$X,
#'                                                                      groups = grouplasso_logreg_data$groups,
#'                                                                      lambda.seq = lambda.seq,
#'                                                                      n.folds = 5,
#'                                                                      b.init.mat = b.mat,
#'                                                                      w = grouplasso_logreg_data$w,
#'                                                                      tol = 1e-3,
#'                                                                      max.iter = 500)
#' @export
grouplasso_logreg_cv_fixedgrid <- function(Y,X,groups,lambda.seq,n.folds,b.init.mat,w,tol=1e-3,max.iter=500)
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
        
        grouplasso_logreg.out <- grouplasso_logreg(rY = Y[-fold.ind],
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
        minus2ll.fold <- - 2 * sum( Y[fold.ind] * log(P.fold) + (1-Y[fold.ind]) * log( 1 - P.fold) )
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


#' Choose tuning parameters by crossvalidation for grouplasso logreg.
#'
#' @param Y the binary response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#'
#' @examples
#' grouplasso_logreg_data <- get_grouplasso_logreg_data(n = 100)
#' 
#' grouplasso_logreg_cv.out <- grouplasso_logreg_cv(Y = grouplasso_logreg_data$Y,
#'                                                  X = grouplasso_logreg_data$X,
#'                                                  groups = grouplasso_logreg_data$groups,
#'                                                  n.lambda = 25,
#'                                                  lambda.min.ratio = 0.01,
#'                                                  n.folds = 5,
#'                                                  w = grouplasso_logreg_data$w,
#'                                                  tol = 1e-3,
#'                                                  max.iter = 500,
#'                                                  report.prog = FALSE)
#' @export
grouplasso_logreg_cv <- function(Y,X,groups,n.lambda,lambda.min.ratio,n.folds,w,tol=1e-4,max.iter=500,report.prog = TRUE){
  
  # obtain lambda.seq from the grid function, as well as the fits on the entire data set,
  # which will be used as initial values for the crossvalidation training fits.
  grouplasso_logreg_grid.out <- grouplasso_logreg_grid(Y = Y,
                                                       X = X,
                                                       groups = groups,
                                                       n.lambda = n.lambda,
                                                       lambda.min.ratio = lambda.min.ratio,
                                                       w = w,
                                                       tol = tol,
                                                       max.iter = max.iter,
                                                       report.prog = report.prog)
                                                               
  lambda.seq <- grouplasso_logreg_grid.out$lambda.seq
  b.mat <- grouplasso_logreg_grid.out$b.mat
  
  # do the crossvalidation
  grouplasso_logreg_cv_fixedgrid.out <- grouplasso_logreg_cv_fixedgrid(Y = Y,
                                                                       X = X,
                                                                       groups = groups,
                                                                       lambda.seq = lambda.seq,
                                                                       n.folds = n.folds,
                                                                       b.init.mat = b.mat,
                                                                       w = w,
                                                                       tol = tol,
                                                                       max.iter = 500)
                                                                               
  output <- list( b.mat = b.mat,
                  b.folds.arr = grouplasso_logreg_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_logreg_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_logreg_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = lambda.seq,
                  iterations = grouplasso_logreg_cv_fixedgrid.out$iterations)
  
  return(output)
  
}


#' Choose tuning parameters by crossvalidation for grouplasso logreg with adaptive weights
#'
#' @param Y the binary response vector
#' @param X matrix containing the design matrices
#' @param groups a vector indicating to which group each covariate belongs
#' @param n.lambda the number of lambda values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @return a list containing the fits over a grid of lambda values as well as the vector of lambda values
#'
#' @examples
#' grouplasso_logreg_data <- get_grouplasso_logreg_data(n = 100)
#' 
#' grouplasso_logreg_cv_adapt.out <- grouplasso_logreg_cv_adapt(Y = grouplasso_logreg_data$Y,
#'                                                              X = grouplasso_logreg_data$X,
#'                                                              groups = grouplasso_logreg_data$groups,
#'                                                              n.lambda = 25,
#'                                                              lambda.min.ratio = 0.01,
#'                                                              n.folds = 5,
#'                                                              w = grouplasso_logreg_data$w,
#'                                                              tol = 1e-3,
#'                                                              max.iter = 500,
#'                                                              report.prog = FALSE)
#' @export
grouplasso_logreg_cv_adapt <- function(Y,X,groups,n.lambda,n.eta,lambda.min.ratio,n.folds,w,tol=1e-3,max.iter=500,report.prog = TRUE){
  
  
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
  grouplasso_logreg.out <- grouplasso_logreg(rY = Y,
                                             rX = X,
                                             groups = groups,
                                             lambda = lambda.initial.fit,
                                             w = w,
                                             tol = tol,
                                             maxiter = max.iter)
                                                     
  
  # now make a new value of w based on this initial fit
  for(j in 1:q){
    
    ind <- which(groups == j)
    w[j] <- min(w[j]/sqrt( sum( grouplasso_logreg.out$beta.hat[ind]^2 )),1e10) # replace Inf with 1e10
    
  }
  
  # obtain lambda.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplasso_logreg_grid.out <- grouplasso_logreg_grid(Y = Y,
                                                       X = X,
                                                       groups = groups,
                                                       n.lambda = n.lambda,
                                                       lambda.min.ratio = lambda.min.ratio,
                                                       w = w,
                                                       tol = tol,
                                                       max.iter = max.iter,
                                                       report.prog = report.prog)
                                                               
  lambda.seq <- grouplasso_logreg_grid.out$lambda.seq
  b.mat <- grouplasso_logreg_grid.out$b.mat
  
  # do the crossvalidation
  grouplasso_logreg_cv_fixedgrid.out <- grouplasso_logreg_cv_fixedgrid(Y = Y,
                                                                       X = X,
                                                                       groups = groups,
                                                                       lambda.seq = lambda.seq,
                                                                       n.folds = n.folds,
                                                                       b.init.mat = b.mat,
                                                                       w = w,
                                                                       tol = tol,
                                                                       max.iter = max.iter)
                                                                               
  output <- list( b.mat = b.mat,
                  b.folds.arr = grouplasso_logreg_cv_fixedgrid.out$b.folds.arr,
                  minus2ll.mat = grouplasso_logreg_cv_fixedgrid.out$minus2ll.mat,
                  which.lambda.cv = grouplasso_logreg_cv_fixedgrid.out$which.lambda.cv,
                  lambda.seq = lambda.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w = w,
                  iterations = grouplasso_logreg_cv_fixedgrid.out$iterations)
  
  return(output)
  
}

