#' Prepare inputs for grouplasso function when using it to fit a semiparametric model
#'
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param d vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @export
semipadd_to_grouplasso <- function(X,nonparm,d,xi,w=1,lambda.beta=1,lambda.f=1){
  
  ww <- ( ( 1 - nonparm ) + nonparm * lambda.f / lambda.beta ) * w
  
  n <- nrow(X)
  pp <- ncol(X)
  
  ##### For first data set
  DD.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length = pp)
  knots.list <- vector("list",length = pp)
  emp.cent <- vector("list",length = pp)
  
  if( length(d) == 1 ){
    
    d <- rep(d,pp)
    
  }
  
  for( j in 1:pp )
  {
    
      
      if(nonparm[j] == 0){
        
        DD.tilde <- cbind(DD.tilde,X[,j])
        groups <- c(groups,j)
        
      } else {
        
        # construct transformed basis function matrices and save important quantities
        
        spsm_cubespline_design.out <- spsm_cubespline_design(X = X[,j],
                                                             d = d[j],
                                                             xi = xi,
                                                             W = NULL)
        
        knots.list[[j]] <- spsm_cubespline_design.out$knots
        emp.cent[[j]] <- spsm_cubespline_design.out$emp.cent
        QQ.inv[[j]] <- spsm_cubespline_design.out$Q.inv
        DD.tilde <- cbind(DD.tilde, spsm_cubespline_design.out$D.tilde )
        
        groups <- c(groups,rep(j,abs(d[j])))
        
      } 
      
  }
  
  output <- list( DD.tilde = DD.tilde,
                  groups = groups,
                  knots.list = knots.list,
                  QQ.inv = QQ.inv,
                  emp.cent = emp.cent,
                  lambda = lambda.beta,
                  w = ww)
  
  return(output)
  
}


#' Convert output from grouplasso to the fitted functions of the semi-parametric additive model
#'
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param groups a vector indicating to which group the entries of the coefficient vector \code{b} belong
#' @param knots.list a list of vectors with the knot locations for nonparametric effects
#' @param emp.cent a list of vectors of the empirical basis function centerings
#' @param QQ.inv the matrix with which to back-transform the group lasso coefficients
#' @param b the group lasso coefficients
#' @export
grouplasso_to_semipadd <- function(X,nonparm,groups,knots.list,emp.cent,QQ.inv,b)
{
  
  n <- nrow(X)
  pp <- ncol(X)
  
  # store fitted functions on data set 1 in a list
  f.hat <- vector("list",pp)
  f.hat.design <- matrix(0,n,pp)
  beta.hat <- rep(NA,pp)
  
  f.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b[1])," }")))
  f.hat.design[,1] <- b[1]
  beta.hat[1] <- b[1]
  
  for(j in 2:pp)
  {
    
    ind <- which(groups == j)
    d <- length(ind)
    
    if(d == 1)
    {
      
      f.hat[[j]] <- eval( parse( text = paste("function(x){ x * ",paste(b[ind])," }")))
      beta.hat[j] <- b[ind]
      
    } else {
      
      Gamma.hat <- QQ.inv[[j]] %*% b[ind]
      f.hat[[j]] <- eval(parse(text = paste("function(x)
                                             {

                                             x <- round(x,10)
                                             x.mat <- spline.des(",paste("c(",paste(round(knots.list[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                             x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
                                             f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma.hat,collapse=","),")",sep=""),")
                                             
                                             return(f.hat)

                                             }"
      )))
      
    }
    
  }
  
  output <- list(f.hat = f.hat,
                 beta.hat = beta.hat)
  
}

#' Fit semiparametric regression model
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
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
#' @examples
#' data <- get_semipadd_data(n = 500, response = "continuous")
#' 
#' semipadd.out <- semipadd(Y = data$Y,
#'                          X = data$X,
#'                          nonparm = data$nonparm,
#'                          response = "continuous",
#'                          w = 1,
#'                          d = 20,
#'                          xi = 1,
#'                          lambda.beta = 1,
#'                          lambda.f = 1,
#'                          tol = 1e-3,
#'                          max.iter = 500)
#' 
#' plot_semipadd(semipadd.out, 
#'               true.functions = list( f = data$f,
#'                                      X = data$X))
#' @export
semipadd <- function(Y,X,nonparm,response,w,d,xi,lambda.beta,lambda.f,tol=1e-4,max.iter=500)
{

  if( any(apply(abs(X),2,sum) == 0)){
    
    stop( "One or more columns of design matrix contain only zeroes")
    
  }
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)
                                            
  if( any(apply(abs(grouplasso_inputs$DD.tilde),2,sum) == 0)){
    
    stop( "A Gram matrix is not positive definite: Try reducing the number of knots." )
    
  }
  
  # get group lasso estimators
  if(response == "continuous"){
    
    grouplasso.out <- grouplasso_linreg(rY = Y,
                                        rX = grouplasso_inputs$DD.tilde,
                                        groups = grouplasso_inputs$groups,
                                        lambda = grouplasso_inputs$lambda,
                                        w = grouplasso_inputs$w,
                                        tol = tol,
                                        maxiter = max.iter)
                                               
  } else if(response == "binary"){
    
    grouplasso.out <- grouplasso_logreg(rY = Y,
                                        rX = grouplasso_inputs$DD.tilde,
                                        groups = grouplasso_inputs$groups,
                                        lambda = grouplasso_inputs$lambda,
                                        w = grouplasso_inputs$w,
                                        tol = tol,
                                        maxiter = max.iter)
                                               
  } else if( response == "gt" ){
    
    grouplasso.out <- grouplasso_gt(Y = Y$I,
                                    Z = Y$A,
                                    Se = Y$Se,
                                    Sp = Y$Sp,
                                    E.approx = Y$E.approx,
                                    X = grouplasso_inputs$DD.tilde,
                                    groups = grouplasso_inputs$groups,
                                    lambda = grouplasso_inputs$lambda,
                                    w = grouplasso_inputs$w,
                                    tol = tol,
                                    maxiter = max.iter)
    
  }
  
  # construct fitted functions from grouplasso output
  semipadd_fitted <- grouplasso_to_semipadd(X = X,
                                            nonparm = nonparm,
                                            groups = grouplasso_inputs$groups,
                                            knots.list = grouplasso_inputs$knots.list,
                                            emp.cent = grouplasso_inputs$emp.cent,
                                            QQ.inv = grouplasso_inputs$QQ.inv,
                                            b = grouplasso.out$beta.hat)
                                                 
  # collect output
  output <- list(f.hat = semipadd_fitted$f.hat,
                 beta.hat = semipadd_fitted$beta.hat,
                 nonparm = nonparm,
                 d = d,
                 xi = xi,
                 knots.list = grouplasso_inputs$knots.list,
                 lambda.beta = lambda.beta,
                 lambda.f = lambda.f)

  class(output) <- "semipadd"

  return(output)

}



#' Compute semiparametric regression model using cv to choose the tuning parameter value
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization for different covariates
#' @param d vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param report.prog a logical indicating whether the progress of the algorithm should be printed to the console
#' @return Returns the estimator of the semiparametric additive model
#' @examples 
#' data <- get_semipadd_data(n = 500, response = "binary")
#' 
#' semipadd_cv_adapt.out <- semipadd_cv_adapt(Y = data$Y,
#'                                            X = data$X,
#'                                            nonparm = data$nonparm,
#'                                            response = "binary",
#'                                            w = 1,
#'                                            d = 20,
#'                                            xi = 1,
#'                                            n.lambda = 5,
#'                                            lambda.min.ratio = .001,
#'                                            lambda.max.ratio = .1,
#'                                            lambda.beta = 1,
#'                                            lambda.f = 1,
#'                                            tol = 1e-3,
#'                                            maxiter = 1000,
#'                                            report.prog = TRUE)
#' 
#' plot_semipadd_cv_adapt(semipadd_cv_adapt.out, 
#'                        true.functions  = list( f = data$f,
#'                                                X = data$X))
#' @export
semipadd_cv_adapt <- function(Y,X,response,nonparm,w,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE){
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)
  
  # get adaptive group lasso estimators under cv choice of lambda
  if( response == "continuous" ){
    
    grouplasso_cv_adapt.out <- grouplasso_linreg_cv_adapt(Y = Y,
                                                          X = grouplasso_inputs$DD.tilde,
                                                          groups = grouplasso_inputs$groups,
                                                          n.lambda = n.lambda,
                                                          lambda.min.ratio = lambda.min.ratio,
                                                          lambda.max.ratio = lambda.max.ratio,
                                                          n.folds = n.folds,
                                                          w = grouplasso_inputs$w,
                                                          tol = tol,
                                                          maxiter = maxiter,
                                                          report.prog = report.prog)
                                                    
  } else if(response == "binary" ){
    
    grouplasso_cv_adapt.out <- grouplasso_logreg_cv_adapt(Y = Y,
                                                          X = grouplasso_inputs$DD.tilde,
                                                          groups = grouplasso_inputs$groups,
                                                          n.lambda = n.lambda,
                                                          lambda.min.ratio = lambda.min.ratio,
                                                          lambda.max.ratio = lambda.max.ratio,
                                                          n.folds = n.folds,
                                                          w = grouplasso_inputs$w,
                                                          tol = tol,
                                                          maxiter = maxiter,
                                                          report.prog = report.prog)
                                                    
  } else if( response == "gt" ){
    
    grouplasso_cv_adapt.out <- grouplasso_gt_cv_adapt(Y = Y$I,
                                                      Z = Y$A,
                                                      Se = Y$Se,
                                                      Sp = Y$Sp,
                                                      E.approx = Y$E.approx,
                                                      X = grouplasso_inputs$DD.tilde,
                                                      groups = grouplasso_inputs$groups,
                                                      n.lambda = n.lambda,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      lambda.max.ratio = lambda.max.ratio,
                                                      n.folds = n.folds,
                                                      w = grouplasso_inputs$w,
                                                      tol = tol,
                                                      maxiter = maxiter,
                                                      report.prog = report.prog)
    
  }
  
  # get matrices of the fitted functions evaluated at the design points
  beta.hat <- matrix(0,ncol(X), n.lambda)
  f.hat <- vector("list",n.lambda)
  
  for(l in 1:n.lambda){
    
    semipadd_fitted <- grouplasso_to_semipadd(X = X,
                                              nonparm = nonparm,
                                              groups = grouplasso_inputs$groups,
                                              knots.list = grouplasso_inputs$knots.list,
                                              emp.cent = grouplasso_inputs$emp.cent,
                                              QQ.inv = grouplasso_inputs$QQ.inv,
                                              b = grouplasso_cv_adapt.out$b.mat[,l])
    
    f.hat[[l]] <- semipadd_fitted$f.hat
    beta.hat[,l] <- semipadd_fitted$beta.hat
    
  }
  
  # prepare output
  output <- list( f.hat = f.hat,
                  beta.hat = beta.hat,
                  nonparm = nonparm,
                  d = d,
                  xi = xi,
                  knots.list = grouplasso_inputs$knots.list,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  n.lambda = n.lambda,
                  n.folds = n.folds,
                  minus2ll.mat = grouplasso_cv_adapt.out$minus2ll.mat,
                  lambda.seq = grouplasso_cv_adapt.out$lambda.seq,
                  lambda.min.ratio = lambda.min.ratio,
                  lambda.max.ratio = lambda.max.ratio,
                  which.lambda.cv = grouplasso_cv_adapt.out$which.lambda.cv,
                  which.lambda.cv.1se = grouplasso_cv_adapt.out$which.lambda.cv.1se,
                  lambda.initial.fit = grouplasso_cv_adapt.out$lambda.initial.fit,
                  iterations = grouplasso_cv_adapt.out$iterations,
                  w = grouplasso_cv_adapt.out$w)
  
  return(output)
  
}


#' Plot method for class semipadd
#' @export
plot_semipadd <- function(x,true.functions = NULL){
  
  f.hat <- x$f.hat
  knots.list <- x$knots.list
  pp <- length(f.hat)
  
  n.plots <- length(which(x$nonparm == 1))
  
  # get evaluations of functions
  f.hat.evals <- matrix(NA,300,pp)
  x.vals <- matrix(NA,300,pp)
  for( j in which(x$nonparm == 1) ){
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    xj.seq <- seq(xj.min,xj.max,length=300)
    x.vals[,j] <- xj.seq
    
    f.hat.evals[,j] <- f.hat[[j]](x.vals[,j])
    
  }
  
  ncols <- 4
  nrows <- ceiling(n.plots/ncols)
  
  par(mfrow=c(nrows,ncols),mar=c(2.1,2.1,1.1,1.1))
  
  for( j in which(x$nonparm == 1) ){
      
    xlims <- c(min(x.vals[,j]),max(x.vals[,j]))
    ylims <- range(f.hat.evals[,-1],na.rm = TRUE)
    
    plot(NA,
         ylim = ylims,
         xlim = xlims)
    
    abline(v=knots.list[[j]],col=rgb(0,0,0,0.15))
    lines(f.hat.evals[,j]~x.vals[,j],col=rgb(0,0,0,1))
    
    if(!is.null(true.functions)){
      
      f.cent.seq <- true.functions$f[[j]](x.vals[,j]) - mean(true.functions$f[[j]](true.functions$X[,j]))
      lines(f.cent.seq ~ x.vals[,j],lty=2)
      
    }
    
  }
    
}

#' Plot method for class semipadd_cv
#' @export
plot_semipadd_cv_adapt <- function(x,true.functions = NULL){
  
  f.hat <- x$f.hat
  knots.list <- x$knots.list
  n.lambda <- x$n.lambda

  nonparm <- x$nonparm
  pp <- length(nonparm)
  
  # get cv choices if they exist
  which.lambda.cv <- x$which.lambda.cv
  which.eta.cv <- x$which.eta.cv
  
  n.plots <- length(which(nonparm == 1))
                    
  # get evaluations of functions
  f.hat.evals <- array(NA,dim=c(300,n.lambda,pp))
  x.vals <- matrix(NA,300,pp)
  for( j in which(x$nonparm == 1) ){
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    xj.seq <- seq(xj.min,xj.max,length = 300)
    x.vals[,j] <- xj.seq
    
    for(l in 1:n.lambda){
  
      f.hat.evals[,l,j] <- f.hat[[l]][[j]](x.vals[,j])
        
    }
    
  }
    
  ncols <- 4
  nrows <- ceiling(n.plots/ncols)
  
  par(mfrow=c(nrows,ncols),mar=c(2.1,2.1,1.1,1.1))
  
  for( j in which(x$nonparm == 1) ){
    
    xlims <- c(min(x.vals[,j]),max(x.vals[,j]))
    ylims <- range(f.hat.evals[,-1,],na.rm = TRUE)
    
    plot(NA,
         ylim = ylims,
         xlim = xlims)
    
    abline(v=knots.list[[j]],col=rgb(0,0,0,0.15))
    
    for(l in 1:n.lambda){
        
      if(length(which.lambda.cv) == 0){
        
        opacity <- 1
          
      } else {
        
        opacity <- ifelse( l == which.lambda.cv,1,0.1)
        
      }
      
      lines(f.hat.evals[,l,j]~x.vals[,j],col=rgb(0,0,0,opacity))
        
    }
    
    if(!is.null(true.functions)){
      
      f.cent.seq <- true.functions$f[[j]](x.vals[,j]) - mean(true.functions$f[[j]](true.functions$X[,j]))
      lines(f.cent.seq ~ x.vals[,j],lty=2)
      
    }
    
  }
                    
}
