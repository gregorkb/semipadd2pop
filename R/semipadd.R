#' Prepare inputs for grouplasso function when using it to fit a semiparametric model
#'
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param d vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param w_int the penalization weights on the interaction terms
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @export
semipadd_to_grouplasso <- function(X,nonparm,d,xi,w=1,int=NULL,w_int=NULL,lambda.beta=1,lambda.f=1)
{
  
  ww <- ( ( 1 - nonparm ) + nonparm * lambda.f / lambda.beta ) * w
  
  n <- nrow(X)
  pp <- ncol(X)
  
  if(length(int) == 0)
  {
    
    ii <- 0
    
  } else {
    
    ii <- nrow(int)
    
  }
  
  ##### For first data set
  DD.tilde <- matrix(NA,n,0)
  groups <- numeric()
  QQ.inv <- vector("list",length = pp + ii)
  knots.list <- vector("list",length = pp + ii )
  emp.cent <- vector("list",length = pp + ii)
  
  if( length(d) == 1 ){
    
    d <- rep(d,pp)
    
  }
  
  for( j in 1:(pp+ii) )
  {
    
    if(j <= pp){ # go through marginal effects
      
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
        
        groups <- c(groups,rep(j,d[j]))
        
      } 
      
    } else if(j > pp){# go through interactions if there are any:
      
        k <- j - pp
      
        if( sum(nonparm[int[k,]]) == 2){
          
          print("Cannot choose two nonparametric components")
          
        } else if(sum(nonparm[int[k,]]) == 1 ){ # int between parametric and nonparametric like f(x)*w
          
          which_nonparm <- int[k,][which(nonparm[int[k,]] == 1)]
          which_parm <- int[k,][which(nonparm[int[k,]] == 0)]
          
          spsm_cubespline_design.out <- spsm_cubespline_design(X = X[,which_nonparm],
                                                               d = d[which_nonparm],
                                                               xi = xi,
                                                               W = X[,which_parm])
          
          knots.list[[j]] <- spsm_cubespline_design.out$knots
          emp.cent[[j]] <- spsm_cubespline_design.out$emp.cent
          QQ.inv[[j]] <- spsm_cubespline_design.out$Q.inv
          DD.tilde <- cbind(DD.tilde, spsm_cubespline_design.out$D.tilde )
          
          groups <- c(groups,rep(j,d[which_nonparm]))
          ww <- c(ww,w_int[k]*lambda.f/lambda.beta)
          
        } else if(sum(nonparm[int[k,]]) == 0){ # int between parametric effects like w1*w2
          
          DD.tilde <- cbind(DD.tilde,X[,int[k,1]]*X[,int[k,2]])
          groups <- c(groups,j)
          ww <- c(ww,w_int[k])
          
        }
        
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
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param groups a vector indicating to which group the entries of the coefficient vector \code{b} belong
#' @param knots.list a list of vectors with the knot locations for nonparametric effects
#' @param emp.cent a list of vectors of the empirical basis function centerings
#' @param QQ.inv the matrix with which to back-transform the group lasso coefficients
#' @param b the group lasso coefficients
#' @export
grouplasso_to_semipadd <- function(X,nonparm,int = NULL,groups,knots.list,emp.cent,QQ.inv,b)
{
  
  n <- nrow(X)
  pp <- ncol(X)
  
  if(length(int) == 0)
  {
    
    ii <- 0
    
  } else {
    
    ii <- nrow(int)
    
  }
  
  # store fitted functions on data set 1 in a list
  f.hat <- vector("list",pp+ii)
  f.hat.design <- matrix(0,n,pp+ii)
  beta.hat <- rep(NA,pp+ii)
  
  f.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b[1])," }")))
  f.hat.design[,1] <- b[1]
  beta.hat[1] <- b[1]
  
  for(j in 2:(pp+ii))
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
    
    if( j <= pp){
      
      f.hat.design[,j] <- f.hat[[j]](X[,j])
    
    # now construct the fitted values for the interaction effects
    } else if(j > pp){
      
      k <- j - pp
      
      if(sum(nonparm[int[k,]]) == 1 ){ # between parametric and nonparametric effect
      
        which_nonparm <- int[k,][which(nonparm[int[k,]] == 1)]
        which_parm <- int[k,][which(nonparm[int[k,]] == 0)]
        
        f.hat.design[,j] <- f.hat[[j]](X[,which_nonparm])*X[,which_parm]
        
      } else if( sum(nonparm[int[k,]]) == 0 ){ # int between parametric effects
        
        f.hat.design[,j] <- b[ind]*X[,int[k,1]]*X[,int[k,2]]
        
      }
      
    }
    
  }
  
  output <- list(f.hat = f.hat,
                 f.hat.design = f.hat.design,
                 beta.hat = beta.hat)
  
}

#' Fit semiparametric regression model
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param w_int the penalization weights on the interaction terms
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
#'                          int = data$int,
#'                          w_int = data$w_int,
#'                          d = 20,
#'                          xi = 1,
#'                          lambda.beta = 1,
#'                          lambda.f = 1,
#'                          tol = 1e-3,
#'                          max.iter = 500)
#' 
#' plot_semipadd(semipadd.out)
#' @export
semipadd <- function(Y,X,nonparm,response,w,int=NULL,w_int=NULL,d,xi,lambda.beta,lambda.f,tol=1e-4,max.iter=500)
{

  if( any(apply(abs(X),2,sum) == 0)){
    
    stop( "One or more columns of design matrix contain only zeroes")
    
  }
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              int = int,
                                              w_int = w_int,
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
                                            int = int,
                                            groups = grouplasso_inputs$groups,
                                            knots.list = grouplasso_inputs$knots.list,
                                            emp.cent = grouplasso_inputs$emp.cent,
                                            QQ.inv = grouplasso_inputs$QQ.inv,
                                            b = grouplasso.out$beta.hat)
                                                 
  # collect output
  output <- list(f.hat = semipadd_fitted$f.hat,
                 f.hat.design = semipadd_fitted$f.hat.design,
                 beta.hat = semipadd_fitted$beta.hat,
                 nonparm = nonparm,
                 d = d,
                 xi = xi,
                 knots.list = grouplasso_inputs$knots.list,
                 lambda.beta = lambda.beta,
                 lambda.f = lambda.f,
                 int = int)

  class(output) <- "semipadd"

  return(output)

}

#' Compute semiparametric regression model while penalizing over a grid of tuning parameter values
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param w covariate-specific weights for different penalization for different covariates
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param w_int the penalization weights on the interaction terms
#' @param d vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' data <- get_semipadd_data(n = 500, response = "binary")
#' 
#' semipadd_grid.out <- semipadd_grid(Y = data$Y,
#'                                    X = data$X,
#'                                    nonparm = data$nonparm,
#'                                    response = "binary",
#'                                    w = 1,
#'                                    int = data$int,
#'                                    w_int = data$w_int,
#'                                    d = 20,
#'                                    xi = 1,
#'                                    n.lambda = 5,
#'                                    lambda.min.ratio = .001,
#'                                    lambda.max.ratio = .1,
#'                                    lambda.beta = 1,
#'                                    lambda.f = 1,
#'                                    tol = 1e-3,
#'                                    maxiter = 500,
#'                                    report.prog = TRUE)
#' 
#' plot_semipadd_grid(semipadd_grid.out)
#' @export
semipadd_grid <- function(Y,X,nonparm,response,w=1,int=NULL,w_int=NULL,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
{
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              int = int,
                                              w_int = w_int,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)
                                                   
                                              
  # get group lasso estimators over a grid of lambda and eta values
  if( response == "continuous"){
    
    grouplasso_grid.out <- grouplasso_linreg_grid(Y = Y,
                                                  X = grouplasso_inputs$DD.tilde,
                                                  groups = grouplasso_inputs$groups,
                                                  n.lambda = n.lambda,
                                                  lambda.min.ratio = lambda.min.ratio,
                                                  lambda.max.ratio = lambda.max.ratio,
                                                  w = grouplasso_inputs$w,
                                                  tol = tol,
                                                  maxiter = maxiter,
                                                  report.prog = report.prog)
                                                         
  
  } else if(response == "binary"){
    
    grouplasso_grid.out <- grouplasso_logreg_grid(Y = Y,
                                                  X = grouplasso_inputs$DD.tilde,
                                                  groups = grouplasso_inputs$groups,
                                                  n.lambda = n.lambda,
                                                  lambda.min.ratio = lambda.min.ratio,
                                                  lambda.max.ratio = lambda.max.ratio,
                                                  w = grouplasso_inputs$w,
                                                  tol = tol,
                                                  maxiter = maxiter,
                                                  report.prog = report.prog)
    
  } else if( response == "gt"){
    
    grouplasso_grid.out <- grouplasso_gt_grid(Y = Y$I,
                                              Z = Y$A,
                                              Se = Y$Se,
                                              Sp = Y$Sp,
                                              E.approx = Y$E.approx,
                                              X = grouplasso_inputs$DD.tilde,
                                              groups = grouplasso_inputs$groups,
                                              n.lambda = n.lambda,
                                              lambda.min.ratio = lambda.min.ratio,
                                              lambda.max.ratio = lambda.max.ratio,
                                              w = grouplasso_inputs$w,
                                              tol = tol,
                                              maxiter = maxiter,
                                              report.prog = report.prog)
                                                  
  }
  
  # get matrices of the fitted functions evaluated at the design points
  
  if(length(int) == 0)
  {
    
    ii <- 0
    
  } else {
    
    ii <- nrow(int)
    
  }
  
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X) + ii,n.lambda))
  beta.hat <- matrix(0,ncol(X) + ii ,n.lambda)
  f.hat <- vector("list",n.lambda)
  
  for(l in 1:n.lambda){
    
    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                int = int,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_grid.out$b.mat[,l])
                                                     
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
                  int = int,
                  lambda.seq = grouplasso_grid.out$lambda.seq,
                  iterations = grouplasso_grid.out$iterations)
  
  class(output) <- "semipadd_grid"
  
  return(output)
  
}


#' Compute semiparametric regression model using cv to choose the tuning parameter value
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param w covariate-specific weights for different penalization for different covariates
#' @param w_int the penalization weights on the interaction terms
#' @param d vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' data <- get_semipadd_data(n = 100, response = "continuous")
#' 
#' semipadd_cv.out <- semipadd_cv(Y = data$Y,
#'                                X = data$X,
#'                                nonparm = data$nonparm,
#'                                response = "continuous",
#'                                w = 1,
#'                                int = data$int,
#'                                w_int = data$w_int,
#'                                d = 20,
#'                                xi = 1,
#'                                n.folds = 5,
#'                                n.lambda = 5,
#'                                lambda.min.ratio = .001,
#'                                lambda.max.ratio = .1,
#'                                lambda.beta = 1,
#'                                lambda.f = 1,
#'                                tol = 1e-3,
#'                                maxiter = 500,
#'                                report.prog = TRUE)
#' 
#' plot_semipadd_grid(semipadd_cv.out)
#' @export
semipadd_cv <- function(Y,X,nonparm,response,int=NULL,w,w_int=NULL,d,xi,n.lambda=5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
{
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              int = int,
                                              w_int = w_int,
                                              lambda.beta = lambda.beta,
                                              lambda.f = lambda.f)
  
  # get group lasso estimators under cv choice of lambda
  if( response == "continuous"){
    
    grouplasso_cv.out <- grouplasso_linreg_cv(Y = Y,
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
                                                  
  } else if(response == "binary"){
    
    grouplasso_cv.out <- grouplasso_logreg_cv(Y = Y,
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
    
  } else if( response == "gt"){
    
    grouplasso_cv.out <- grouplasso_gt_cv(Y = Y$I,
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
  if(length(int) == 0)
  {
    
    ii <- 0
    
  } else {
    
    ii <- nrow(int)
    
  }
  
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X) + ii, n.lambda))
  beta.hat <- matrix(0,ncol(X) + ii, n.lambda)
  f.hat <- vector("list",n.lambda)
  
  for(l in 1:n.lambda){
    
    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                int = int,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_cv.out$b.mat[,l])
    
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
                  n.folds = n.folds,
                  lambda.min.ratio = lambda.min.ratio,
                  lambda.max.ratio = lambda.max.ratio,
                  lambda.seq = grouplasso_cv.out$lambda.seq,
                  which.lambda.cv = grouplasso_cv.out$which.lambda.cv,
                  iterations = grouplasso_cv.out$iterations,
                  int = int)
  
  class(output) <- "sempadd_cv"
  
  return(output)
  
}


#' Compute semiparametric regression model using cv to choose the tuning parameter value
#'
#' @param Y the response data
#' @param X the matrix with the observed covariate values (including a column of ones for the intercept)
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param nonparm a vector indicating for which covariates a nonparametric function is to be estimated
#' @param int a matrix with rows giving the pairs of covariates for which to include interaction terms
#' @param w covariate-specific weights for different penalization for different covariates
#' @param w_int the penalization weights on the interaction terms
#' @param d vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param n.lambda the number of lambda values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#' @examples 
#' data <- get_semipadd_data(n = 500, response = "binary")
#' 
#' semipadd_cv_adapt.out <- semipadd_cv_adapt(Y = data$Y,
#'                                            X = data$X,
#'                                            nonparm = data$nonparm,
#'                                            response = "binary",
#'                                            w = 1,
#'                                            int = data$int,
#'                                            w_int = data$w_int,
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
#' plot_semipadd_grid(semipadd_cv_adapt.out)
#' @export
semipadd_cv_adapt <- function(Y,X,response,nonparm,int=NULL,w,w_int=NULL,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
{
  
  # prepare input for grouplasso function
  grouplasso_inputs <- semipadd_to_grouplasso(X = X,
                                              nonparm = nonparm,
                                              d = d,
                                              xi = xi,
                                              w = w,
                                              int = int,
                                              w_int = w_int,
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
  if(length(int) == 0)
  {
    
    ii <- 0
    
  } else {
    
    ii <- nrow(int)
    
  }
  
  f.hat.design <- array(0,dim=c(nrow(X),ncol(X) + ii, n.lambda))
  beta.hat <- matrix(0,ncol(X) + ii, n.lambda)
  f.hat <- vector("list",n.lambda)
  
  for(l in 1:n.lambda){
    
    semipadd_fitted <- grouplasso_to_semipadd(X = X,
                                              nonparm = nonparm,
                                              int = int,
                                              groups = grouplasso_inputs$groups,
                                              knots.list = grouplasso_inputs$knots.list,
                                              emp.cent = grouplasso_inputs$emp.cent,
                                              QQ.inv = grouplasso_inputs$QQ.inv,
                                              b = grouplasso_cv_adapt.out$b.mat[,l])
    
    f.hat[[l]] <- semipadd_fitted$f.hat
    f.hat.design[,,l] <- semipadd_fitted$f.hat.design
    beta.hat[,l] <- semipadd_fitted$beta.hat
    
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
                  n.folds = n.folds,
                  lambda.seq = grouplasso_cv_adapt.out$lambda.seq,
                  lambda.min.ratio = lambda.min.ratio,
                  lambda.max.ratio = lambda.max.ratio,
                  which.lambda.cv = grouplasso_cv_adapt.out$which.lambda.cv,
                  lambda.initial.fit = grouplasso_cv_adapt.out$lambda.initial.fit,
                  iterations = grouplasso_cv_adapt.out$iterations,
                  int = int)
  
  class(output) <- "sempadd_cv"
  
  return(output)
  
}



    