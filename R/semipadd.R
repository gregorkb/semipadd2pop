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
        
        if(d[j] < 0){
          
          int.knots <- seq(min(X[,j]), max(X[,j]), length = -d[j] - 2 + 1 )
          d[j] <- -d[j]
          
        } else {
          
          int.knots <- quantile(X[,j],seq(0,1,length=d[j]-2+1)) # add one, so that one can be removed after centering to restore full-rank.
          
        }
        
        boundary.knots <- range(int.knots)
        all.knots <- sort(c(rep(boundary.knots,3),int.knots))
        knots.list[[j]] <- all.knots
        
        Bj <- spline.des(all.knots,X[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
        emp.cent[[j]] <- apply(Bj,2,mean)
        Bj.cent <- Bj - matrix(emp.cent[[j]],n,d[j],byrow=TRUE)
        
        # construct matrix in which l2 norm of function is a quadratic form
        M <- t(Bj.cent) %*% Bj.cent / n
        
        # construct matrix in which 2nd derivative penalty is a quadratic form
        R <- matrix(NA,d[j]+1,d[j]+1)
        dsq_bspline.mat <- spline.des(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)$design
        for(k in 1:(d[j]+1))
          for(l in 1:(d[j]+1))
          {
            
            pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
            h <- diff(int.knots)
            R[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.
            
          }
        
        Q <- chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
        Q.inv <- solve(Q)
        
        QQ.inv[[j]] <- Q.inv
        
        # construct transformed basis function matrices
        DD.tilde <- cbind(DD.tilde, Bj.cent %*% Q.inv)
        
        groups <- c(groups,rep(j,d[j]))
        
      } 
      
    } else if(j > pp){# go through interactions if there are any:
      
        k <- j - pp
      
        if( sum(nonparm[int[k,]]) == 2){
          
          print("Cannot choose two nonparametric components")
          
        } else if(sum(nonparm[int[k,]]) == 1 ){ # int between parametric and nonparametric like f(x)*w
          
          which_nonparm <- int[k,][which(nonparm[int[k,]] == 1)]
          which_parm <- int[k,][which(nonparm[int[k,]] == 0)]
          
          if(d[which_nonparm] < 0){
            
            int.knots <- seq(min(X[,which_nonparm]), max(X[,which_nonparm]), length = -d[which_nonparm] - 2 + 1 )
            d[which_nonparm] <- -d[which_nonparm]
            
          } else {
            
            int.knots <- quantile(X[,which_nonparm],seq(0,1,length=d[which_nonparm]-2+1)) # add one, so that one can be removed after centering to restore full-rank.
            
          }
          
          boundary.knots <- range(int.knots)
          all.knots <- sort(c(rep(boundary.knots,3),int.knots))
          knots.list[[j]] <- all.knots
          
          # build interaction design matrix
          Bj <- diag(X[,which_parm]) %*% spline.des(all.knots,X[,which_nonparm],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
          
          emp.cent.int <- apply(Bj[which(X[,which_parm]==1),],2,mean)
          Bj.cent <- Bj - matrix(emp.cent.int,n,d[which_nonparm],byrow=TRUE)
          
          # construct matrix in which l2 norm of function is a quadratic form
          M <- t(Bj.cent) %*% Bj.cent / n
          
          # construct matrix in which 2nd derivative penalty is a quadratic form
          R <- matrix(NA,d[which_nonparm]+1,d[which_nonparm]+1)
          dsq_bspline.mat <- spline.des(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)$design
          for(m in 1:(d[which_nonparm]+1))
            for(l in 1:(d[which_nonparm]+1))
            {
              
              pcwiselin <- dsq_bspline.mat[,m] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
              h <- diff(int.knots)
              R[m,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.
              
            }
          
          Q <-  chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
          Q.inv <- solve(Q)
          
          # construct transformed basis function matrices
          DD.tilde <- cbind(DD.tilde, Bj.cent %*% Q.inv)
          QQ.inv[[j]] <- Q.inv
          emp.cent[[j]] <- emp.cent.int
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
  
  f.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b[1])," }")))
  f.hat.design[,1] <- b[1]
  
  beta.hat <- rep(NA,pp)
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
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt}
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
#'
#' @examples
#' data <- get_semipadd_linreg_wint_data(n = 2000)
#' 
#' semipadd_linreg_wint.out <- semipadd_linreg_wint(Y = data$Y,
#'                                                  X = data$X,
#'                                                  nonparm = data$nonparm,
#'                                                  w = 1,
#'                                                  int = data$int,
#'                                                  w_int = data$w_int,
#'                                                  d = 20,
#'                                                  xi = 1,
#'                                                  lambda.beta = 1,
#'                                                  lambda.f = 1,
#'                                                  tol = 1e-4,
#'                                                  max.iter = 500)
#' 
#' plot_semipaddgt_wint(semipadd_linreg_wint.out)
#' 
#' semipadd_linreg_wint.out$beta.hat
#' @export
semipadd <- function(Y,X,nonparm,response,w,int=NULL,w_int=NULL,d,xi,lambda.beta,lambda.f,tol=1e-4,max.iter=500)
{

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
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt}
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
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt}
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
#' @export
semipadd_cv <- function(Y,X,nonparm,response,int,w,w_int,d,xi,n.lambda=5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt}
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
#' @export
semipadd_cv_adapt <- function(Y,X,response,nonparm,int,w,w_int,d,xi,n.lambda = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
    
    semipaddgt_fitted <- grouplasso_to_semipadd(X = X,
                                                nonparm = nonparm,
                                                int = int,
                                                groups = grouplasso_inputs$groups,
                                                knots.list = grouplasso_inputs$knots.list,
                                                emp.cent = grouplasso_inputs$emp.cent,
                                                QQ.inv = grouplasso_inputs$QQ.inv,
                                                b = grouplasso_cv_adapt.out$b.mat[,l])
    
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
                  lambda.seq = grouplasso_cv_adapt.out$lambda.seq,
                  lambda.min.ratio = lambda.min.ratio,
                  lambda.max.ratio = lambda.max.ratio,
                  which.lambda.cv = grouplasso_cv_adapt.out$which.lambda.cv,
                  lambda.initial.fit = grouplasso_cv_adapt.out$lambda.initial.fit,
                  iterations = grouplasso_cv_adapt.out$iterations)
  
  class(output) <- "sempadd_cv"
  
  return(output)
  
}



#' Plot method for class semipadd
#' @export
plot_semipadd <- function(x)
{
  
  f.hat <- x$f.hat
  f.hat.design <- x$f.hat.design
  knots.list <- x$knots.list
  nonparm <- x$nonparm
  pp <- length(nonparm)
  n.plots <- length(which(nonparm == 1))
  int <- x$int
  
  ncols <- 4
  nrows <- ceiling(n.plots/ncols)
  
  par(mfrow = c(nrows,ncols), mar = c(2.1,2.1,1.1,1.1))
  
  for( j in which(nonparm == 1) ){
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    
    xj.seq <- seq(xj.min,xj.max,length=200)
    
    plot(NA,ylim = range(f.hat.design[,-1]),xlim=c(xj.min,xj.max))
    if(nonparm[j]==1) abline(v=knots.list[[j]],col=rgb(0,0,0,0.15))
    
    
    
    plot(f.hat[[j]],xj.min,xj.max,add=TRUE,col=rgb(0,0,0,1))
    
    if(length(x$int)!=0){
      
      if(any(int == j)){
        
        which.interactions <- which(int == j, arr.ind = TRUE)[,1]
        
        for( k in (which.interactions + pp))
        {
          y.seq <- f.hat[[k]](xj.seq) + f.hat[[j]](xj.seq)
          lines(y.seq~xj.seq,col=rgb(0,0,0,1))
          
        }
      }
    }
  }
}

#' Plot method for class semipadd_grid
#' @export
plot_semipadd_grid <- function(x)
{
  
  f.hat <- x$f.hat
  f.hat.design <- x$f.hat.design
  knots.list <- x$knots.list
  nonparm <- x$nonparm
  pp <- length(nonparm)
  n.plots <- length(which(nonparm == 1))
  int <- x$int
  n.lambda <- x$n.lambda
  
  ncols <- 4
  nrows <- ceiling(n.plots/ncols)
  
  par(mfrow = c(nrows,ncols), mar = c(2.1,2.1,1.1,1.1))
  
  for( j in which(nonparm == 1) ){
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    
    xj.seq <- seq(xj.min,xj.max,length=200)
    
    plot(NA,ylim = range(f.hat.design[,-1,]),xlim=c(xj.min,xj.max))
    if(nonparm[j]==1) abline(v=knots.list[[j]],col=rgb(0,0,0,0.15))
    
    
    for(l in 1:n.lambda){
      
      plot(f.hat[[l]][[j]],xj.min,xj.max,add=TRUE,col=rgb(0,0,0,.5))
      
      if(length(x$int)!=0){
        
        if(any(int == j)){
          
          which.interactions <- which(int == j, arr.ind = TRUE)[,1]
          
          for( k in (which.interactions + pp))
          {
            y.seq <- f.hat[[l]][[k]](xj.seq) + f.hat[[l]][[j]](xj.seq)
            lines(y.seq~xj.seq,col=rgb(0,0,0,.5))
            
          }
        }
      }
    }
  }
}


#' Plot method for class semipadd_cv
#' @export
plot_semipadd_cv <- function(x)
{
  
  f.hat <- x$f.hat
  f.hat.design <- x$f.hat.design
  knots.list <- x$knots.list
  nonparm <- x$nonparm
  pp <- length(nonparm)
  n.plots <- length(which(nonparm == 1))
  int <- x$int
  n.lambda <- x$n.lambda
  which.lambda.cv <- x$which.lambda.cv
  
  ncols <- 4
  nrows <- ceiling(n.plots/ncols)
  
  par(mfrow = c(nrows,ncols), mar = c(2.1,2.1,1.1,1.1))
  
  for( j in which(nonparm == 1) ){
    
    xj.min <- min(knots.list[[j]]) + 1e-2
    xj.max <- max(knots.list[[j]]) - 1e-2
    
    xj.seq <- seq(xj.min,xj.max,length=200)
    
    plot(NA,ylim = range(f.hat.design[,-1,]),xlim=c(xj.min,xj.max))
    if(nonparm[j]==1) abline(v=knots.list[[j]],col=rgb(0,0,0,0.15))
    
    for(l in 1:n.lambda){
      
      opacity <- ifelse(l == which.lambda.cv,1,0.25)
      plot(f.hat[[l]][[j]],xj.min,xj.max,add=TRUE,col=rgb(0,0,0,opacity))
      
      if(length(x$int)!=0){
        
        if(any(int == j)){
          
          which.interactions <- which(int == j, arr.ind = TRUE)[,1]
          
          for( k in (which.interactions + pp))
          {
            y.seq <- f.hat[[l]][[k]](xj.seq) + f.hat[[l]][[j]](xj.seq)
            lines(y.seq~xj.seq,col=rgb(0,0,0,opacity))
            
          }
        }
      }
    }
  }
}
    



#' Generate a data set for semiparametric additive modeling with continuous responses that has an interaction effect
#'
#' @param n the sample size
#' @return a list containing the data
#' @export
get_semipadd_data <- function(n,response = "continuous")
{
  
  p <- 6
  q <- 4
  zeta1 <- 3/20
  zeta2 <- 10/20
  W <- cbind(corrBern(n,probs=c(1:p)/(2*p),Rho = zeta1^abs(outer(1:p,1:p,"-"))))
  X <- (corrUnif(n,Rho = zeta2^abs(outer(1:q,1:q,"-")))-.5)*5
  
  XX <- cbind(1,W[,c(1,2)],X[,c(1,2)],W[,-c(1,2)],X[,-c(1,2)])
  nonparm <- c(0,rep(0,2),rep(1,2),rep(0,p - 2),rep(1,q - 2))
  pp <- ncol(XX)
  
  # set up the true functions
  f <- vector("list",11)
  f[[1]] <- function(x){0} # intercept
  f[[2]] <- function(x){x*2}
  f[[3]] <- function(x){x*0}
  f[[4]] <- function(x){-2*sin(x*2)}
  f[[5]] <- function(x){-x}
  f[[6]] <- function(x){x*2}
  f[[7]] <- function(x){x*0}
  f[[8]] <- function(x){0*x}
  f[[9]] <- function(x){0*x}
  f[[10]] <- function(x){0*x}
  f[[11]] <- function(x){(exp(-x)-2/5*sinh(5/2))/2}
  
  # record coefficients for covariates to be fit parametrically
  beta <- c(0,2,0,NA,NA,2,0,0,0,NA,NA)
  
  # make marginal effects
  f.design <- matrix(0,n,11)
  for(j in 1:11){
    
    if(nonparm[j] == 1){
      
      f.design[,j] <- f[[j]](XX[,j]) - mean(f[[j]](XX[,j]))
      
    } else {
      
      f.design[,j] <- f[[j]](XX[,j])
      
    }
    
  }
  
  # make nonparametric interaction effect:
  f_int <- function(x){2*(pnorm(x - 1) - .5)}
  nonparm_int_effect <- (f_int(XX[,4]) - mean(f_int(XX[which(XX[,6] == 1),4])))*XX[,6]
  
  # make parametric interaction effect:
  
  parm_int_effect <- - 1 * XX[,2] * XX[,6]
  
  if( response == "continuous")
  {
    
    Y <- apply(f.design,1,sum) + nonparm_int_effect + parm_int_effect + rnorm(n,0,.5)
    
  } else if(response == "binary"){
    
    P <- logit(apply(f.design,1,sum) + nonparm_int_effect + parm_int_effect)
    Y <- rbinom(n,1,prob=P)
    
  } else if( response == "gt"){
    
    P <- logit(apply(f.design,1,sum) + nonparm_int_effect + parm_int_effect)
    Y.true <- rbinom(n,1,prob=P)
    
    Se <- c(.98,.96)
    Sp <- c(.97,.99)
    assay.out <- dorfman.assay.gen(Y.true,Se,Sp,cj=4)
    
    Y <- list(  A = assay.out$Z,
                I = assay.out$Y,
                Se = Se,
                Sp = Sp,
                E.approx = FALSE)
    
  } else {
    
    stop("invalid response type")
    
  }
  
  
  int = matrix(c(4,6,
                 4,7,
                 5,6,
                 2,6),
               nrow = 4, 
               byrow=TRUE)
  
  w_int = c(1,1,1,1)
  
  output <- list(X = XX,
                 nonparm = nonparm,
                 f = f,
                 beta = beta,
                 Y = Y,
                 int = int,
                 w_int = w_int)
  
  return(output)
  
}
    