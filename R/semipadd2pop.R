
#' Prepare inputs for grouplasso2pop function when using it to fit a semiparametric model
#'
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}.
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @export
semipadd2pop_to_grouplasso2pop <- function(X1,nonparm1,X2,nonparm2,nCom,d1,d2,xi,w1=1,w2=1,w=1,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1)
{
  
  ww1 <- ((1-nonparm1) + nonparm1 * lambda.f/lambda.beta) * w1
  ww2 <- ((1-nonparm2) + nonparm2 * lambda.f/lambda.beta) * w2
  ww1[1] <- 0 # do not penalize intercept
  ww2[1] <- 0 # do not penalize intercept
  ww <- ((1-nonparm1) + nonparm1 * eta.f/eta.beta) * w
  Com <- 1 + 1:nCom
  ww <- ww[1:max(Com)]
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  pp1 <- ncol(X1)
  pp2 <- ncol(X2)
  
  if( length(d1) == 1 ){
    
    d1 <- rep(d1,pp1)
    
  }
  
  if( length(d2) == 1 ){
    
    d2 <- rep(d2,pp1)
    
  }
  
  ##### For first data set
  DD1.tilde <- matrix(NA,n1,0)
  groups1 <- numeric()
  QQ1.inv <- vector("list",length=pp1)
  knots.list1 <- vector("list",length=pp1)
  emp.cent1 <- vector("list",length=pp1)
  
  for( j in 1:pp1 ){
    
    if(nonparm1[j] == 0){
      
      DD1.tilde <- cbind(DD1.tilde,X1[,j])
      groups1 <- c(groups1,j)
      
    } else {
      
      spsm_cubespline_design.out <- spsm_cubespline_design(X = X1[,j],
                                                           d = d1[j],
                                                           xi = xi,
                                                           W = NULL)
      
      knots.list1[[j]] <- spsm_cubespline_design.out$knots
      emp.cent1[[j]] <- spsm_cubespline_design.out$emp.cent
      QQ1.inv[[j]] <- spsm_cubespline_design.out$Q.inv
      DD1.tilde <- cbind(DD1.tilde, spsm_cubespline_design.out$D.tilde)
      d1[j] <- spsm_cubespline_design.out$d
      
      groups1 <- c(groups1,rep(j,d1[j]))
      
    }
    
  }
  
  ##### For second data set
  DD2.tilde <- matrix(NA,n2,0)
  groups2 <- numeric()
  QQ2.inv <- vector("list",length=pp2)
  knots.list2 <- vector("list",length=pp2)
  emp.cent2 <- vector("list",length=pp2)
  
  for( j in 1:pp2 )
  {
    
    if(nonparm2[j] == 0){
      
      DD2.tilde <- cbind(DD2.tilde,X2[,j])
      groups2 <- c(groups2,j)
      
    } else {
      
      spsm_cubespline_design.out <- spsm_cubespline_design(X = X2[,j],
                                                           d = d2[j],
                                                           xi = xi,
                                                           W = NULL)
      
      knots.list2[[j]] <- spsm_cubespline_design.out$knots
      emp.cent2[[j]] <- spsm_cubespline_design.out$emp.cent
      QQ2.inv[[j]] <- spsm_cubespline_design.out$Q.inv
      DD2.tilde <- cbind(DD2.tilde, spsm_cubespline_design.out$D.tilde)
      d2[j] <- spsm_cubespline_design.out$d
        
      groups2 <- c(groups2,rep(j,d2[j]))
      
    }
    
  }
  
  
  # now construct matrices needed for the dissimilarity penalties for the effects of common covariates
  AA1.tilde <- vector("list",length=pp1)
  AA2.tilde <- vector("list",length=pp2)
  
  if(length(Com)>0){
    for( j in Com ){
      
      if(nonparm1[j]==0){ # nonparm1 and nonparm2 are identical over the indices in Com
        
        AA1.tilde[[j]] <- as.matrix(1)
        AA2.tilde[[j]] <- as.matrix(1)
        
      } else if(nonparm1[j]==1){
        
        X1j.srt <- sort(X1[,j])
        X2j.srt <- sort(X2[,j])
        
        Xj.union <- sort(c(X1[,j],X2[,j]))
        int.ind <- which( max(X1j.srt[1],X2j.srt[1]) < Xj.union & Xj.union < min(X1j.srt[n1],X2j.srt[n2]) )
        
        Xj.int <- Xj.union[int.ind]
        
        AA1j <- spline.des(knots.list1[[j]],Xj.int,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove same one
        AA1j.cent <- AA1j - matrix(apply(AA1j,2,mean),length(int.ind),d1[j],byrow=TRUE)
        AA1.tilde[[j]] <- AA1j.cent %*% QQ1.inv[[j]]
        
        AA2j <- spline.des(knots.list2[[j]],Xj.int,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove same one
        AA2j.cent <- AA2j - matrix(apply(AA2j,2,mean),length(int.ind),d2[j],byrow=TRUE)
        AA2.tilde[[j]] <- AA2j.cent %*% QQ2.inv[[j]]
        
      }
    }
  }
  

  output <- list( DD1.tilde = DD1.tilde,
                  DD2.tilde = DD2.tilde,
                  groups1 = groups1,
                  groups2 = groups2,
                  knots.list1 = knots.list1,
                  knots.list2 = knots.list2,
                  QQ1.inv = QQ1.inv,
                  QQ2.inv = QQ2.inv,
                  emp.cent1 = emp.cent1,
                  emp.cent2 = emp.cent2,
                  AA1.tilde = AA1.tilde,
                  AA2.tilde = AA2.tilde,
                  lambda = lambda.beta,
                  eta = eta.beta,
                  w1 = ww1,
                  w2 = ww2,
                  w = ww,
                  Com = Com)
  
  return(output)
  
}

#' Convert output from grouplasso2pop to the fitted functions of the semi-parametric additive model
#'
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param groups1 a vector indicating to which group the entries of the coefficient vector \code{b1} belong
#' @param knots.list1 a list of vectors with the knot locations for nonparametric effects for data set 1
#' @param emp.cent1 a list of vectors of the empirical basis function centerings for data set 1
#' @param QQ1.inv the matrix with which to back-transform the group lasso coefficients for data set 1
#' @param b1 the group lasso coefficients for data set 1
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param groups2 a vector indicating to which group the entries of the coefficient vector \code{b2} belong
#' @param knots.list2 a list of vectors with the knot locations for nonparametric effects for data set 2
#' @param emp.cent2 a list of vectors of the empirical basis function centerings for data set 2
#' @param QQ2.inv the matrix with which to back-transform the group lasso coefficients for data set 2
#' @param b2 the group lasso coefficients for data set 2
#' @export
grouplasso2pop_to_semipadd2pop <- function(X1,nonparm1,groups1,knots.list1,emp.cent1,QQ1.inv,b1,X2,nonparm2,groups2,knots.list2,emp.cent2,QQ2.inv,b2)
{
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  pp1 <- ncol(X1)
  pp2 <- ncol(X2)
  
  # store fitted functions on data set 1 in a list
  f1.hat <- vector("list",pp1)
  f1.hat.design <- matrix(0,n1,pp1)
  
  f1.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b1[1])," }")))
  f1.hat.design[,1] <- b1[1]
  
  
  beta1.hat <- rep(NA,pp1)
  beta1.hat[1] <- b1[1]
  for(j in 2:pp1)
  {
    
    if(nonparm1[j] == 0)
    {
      
      ind <- which(groups1 == j)
      f1.hat[[j]] <- eval( parse( text = paste("function(x){ x * ",paste(b1[ind])," }")))
      beta1.hat[j] <- b1[ind]
      
    } else {
      
      ind <- which(groups1 == j)
      d <- length(ind)
      
      Gamma1.hat <- QQ1.inv[[j]] %*% b1[ind]
      f1.hat[[j]] <- eval(parse(text = paste("function(x)
                                             {

                                             x <- round(x,10)
                                             x.mat <- spline.des(",paste("c(",paste(round(knots.list1[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                             x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent1[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
                                             f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma1.hat,collapse=","),")",sep=""),")

                                             return(f.hat)

                                             }"
      )))
      
    }
    
    f1.hat.design[,j] <- f1.hat[[j]](X1[,j])
    
  }
  
  # store fitted functions on data set 2 in a list
  f2.hat <- vector("list",pp2)
  f2.hat.design <- matrix(0,n2,pp2)
  
  f2.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b2[1])," }")))
  f2.hat.design[,1] <- b2[1]
  
  beta2.hat <- rep(NA,pp2)
  beta2.hat[1] <- b2[1]
  for(j in 2:pp2)
  {
    
    if(nonparm2[j] == 0)
    {
      
      ind <- which(groups2 == j)
      f2.hat[[j]] <- eval( parse( text = paste("function(x){ x * ",paste(b2[ind])," }")))
      beta2.hat[j] <- b2[ind]
      
    } else {
      
      ind <- which(groups2 == j)
      d <- length(ind)
      
      Gamma2.hat <- QQ2.inv[[j]] %*% b2[ind]
      f2.hat[[j]] <- eval(parse(text=paste("function(x)
                                           {

                                           x <- round(x,10)
                                           x.mat <- spline.des(",paste("c(",paste(round(knots.list2[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                           x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent2[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
                                           f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma2.hat,collapse=","),")",sep=""),")

                                           return(f.hat)

                                           }"
      )))
      
    }
    
    f2.hat.design[,j] <- f2.hat[[j]](X2[,j])
    
  }
  
  output <- list(f1.hat = f1.hat,
                 f1.hat.design = f1.hat.design,
                 f2.hat = f2.hat,
                 f2.hat.design = f2.hat.design,
                 beta1.hat = beta1.hat,
                 beta2.hat = beta2.hat)
  
}


#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 the response data from data set 1
#' @param XX1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the response data from data set 2
#' @param XX2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
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
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' data <- get_semipadd2pop_data(n1 = 500, n2 = 600, model = 1, response = "gt")
#' 
#' semipadd.out <- semipadd2pop(Y1 = data$Y1,
#'                              X1 = data$X1,
#'                              nonparm1 = data$nonparm1,
#'                              Y2 = data$Y2,
#'                              X2 = data$X2,
#'                              nonparm2 = data$nonparm2,
#'                              response = "gt",
#'                              rho1 = 2,
#'                              rho2 = 1,
#'                              nCom = data$nCom,
#'                              d1 = 25,
#'                              d2 = 15,
#'                              xi = .5,
#'                              w1 = 1,
#'                              w2 = 1,
#'                              w = 1,
#'                              lambda.beta = .01,
#'                              lambda.f = .01,
#'                              eta.beta = .1, 
#'                              eta.f = .1,
#'                              tol = 1e-3,
#'                              maxiter = 500)
#' 
#' plot_semipadd2pop_gt(semipadd.out)
#' @export
semipadd2pop <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,response,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,tol=1e-4,maxiter=500,plot_obj=FALSE)
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
                                                          w1 = w1,
                                                          w2 = w2,
                                                          w = w,
                                                          lambda.beta = lambda.beta,
                                                          lambda.f = lambda.f,
                                                          eta.beta = eta.beta,
                                                          eta.f = eta.f)
  
  # get group lasso estimators
  if( response == "continuous"){
    
    grouplasso2pop.out <- grouplasso2pop_linreg(rY1 = Y1,
                                                rX1 = grouplasso2pop_inputs$DD1.tilde,
                                                groups1 = grouplasso2pop_inputs$groups1,
                                                rY2 = Y2,
                                                rX2 = grouplasso2pop_inputs$DD2.tilde,
                                                groups2 = grouplasso2pop_inputs$groups2,
                                                rho1 = rho1,
                                                rho2 = rho2,
                                                lambda = grouplasso2pop_inputs$lambda,
                                                eta = grouplasso2pop_inputs$eta,
                                                w1 = grouplasso2pop_inputs$w1,
                                                w2 = grouplasso2pop_inputs$w2,
                                                w = grouplasso2pop_inputs$w,
                                                rAA1 = grouplasso2pop_inputs$AA1.tilde,
                                                rAA2 = grouplasso2pop_inputs$AA2.tilde,
                                                rCom = grouplasso2pop_inputs$Com,
                                                tol = tol,
                                                maxiter = maxiter)
                                                       
  
  } else if( response == "binary" ){
    
    grouplasso2pop.out <- grouplasso2pop_logreg(rY1 = Y1,
                                                rX1 = grouplasso2pop_inputs$DD1.tilde,
                                                groups1 = grouplasso2pop_inputs$groups1,
                                                rY2 = Y2,
                                                rX2 = grouplasso2pop_inputs$DD2.tilde,
                                                groups2 = grouplasso2pop_inputs$groups2,
                                                rho1 = rho1,
                                                rho2 = rho2,
                                                lambda = grouplasso2pop_inputs$lambda,
                                                eta = grouplasso2pop_inputs$eta,
                                                w1 = grouplasso2pop_inputs$w1,
                                                w2 = grouplasso2pop_inputs$w2,
                                                w = grouplasso2pop_inputs$w,
                                                rAA1 = grouplasso2pop_inputs$AA1.tilde,
                                                rAA2 = grouplasso2pop_inputs$AA2.tilde,
                                                rCom = grouplasso2pop_inputs$Com,
                                                tol = tol,
                                                maxiter = maxiter)
                                                       
  } else if(response == "gt"){
    
    grouplasso2pop.out <- grouplasso2pop_gt(Y1 = Y1$I,
                                            Z1 = Y1$A,
                                            Se1 = Y1$Se,
                                            Sp1 = Y1$Sp,
                                            E.approx1 = Y1$E.approx,
                                            X1 = grouplasso2pop_inputs$DD1.tilde,
                                            groups1 = grouplasso2pop_inputs$groups1,
                                            Y2 = Y2$I,
                                            Z2 = Y2$A,
                                            Se2 = Y2$Se,
                                            Sp2 = Y2$Sp,
                                            E.approx2 = Y2$E.approx,
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
                                            tol = tol,
                                            maxiter = maxiter)
                                                       
  }
  
  # construct fitted functions from grouplasso2pop output
  semipadd2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                        nonparm1 = nonparm1,
                                                        groups1 = grouplasso2pop_inputs$groups1,
                                                        knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                        emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                        QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                        b1 = grouplasso2pop.out$beta1.hat,
                                                        X2 = X2,
                                                        nonparm2 = nonparm2,
                                                        groups2 = grouplasso2pop_inputs$groups2,
                                                        knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                        emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                        QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                        b2 = grouplasso2pop.out$beta2.hat)
  
  # collect output
  output <- list( f1.hat = semipadd2pop_fitted$f1.hat,
                  f2.hat = semipadd2pop_fitted$f2.hat,
                  f1.hat.design = semipadd2pop_fitted$f1.hat.design,
                  f2.hat.design = semipadd2pop_fitted$f2.hat.design,
                  beta1.hat = semipadd2pop_fitted$beta1.hat,
                  beta2.hat = semipadd2pop_fitted$beta2.hat,
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
                  Com = grouplasso2pop_inputs$Com)
  
  class(output) <- "semipadd2pop"
  
  return(output)
}

#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity over a grid of tuning parameter values
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the binary response vector of data set 2
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param n.lambda the number of lambda values with which to make the grid
#' @param n.eta the number of eta values with which to make the grid
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param eta.min.ratio ratio of the smallest to largest value in the sequence of eta values
#' @param eta.max.ratio controls the largest value of eta in the eta sequence
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#'                                                          n2 = 604)
#'
#' semipadd2pop_logreg_grid.out <- semipadd2pop_logreg_grid(Y1 = semipadd2pop_logreg_data$Y1,
#'                                                          X1 = semipadd2pop_logreg_data$X1,
#'                                                          nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                          Y2 = semipadd2pop_logreg_data$Y2,
#'                                                          X2 = semipadd2pop_logreg_data$X2,
#'                                                          nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                          rho1 = 2,
#'                                                          rho2 = 1,
#'                                                          nCom = semipadd2pop_logreg_data$nCom,
#'                                                          d1 = semipadd2pop_logreg_data$nonparm1*25,
#'                                                          d2 = semipadd2pop_logreg_data$nonparm2*15,
#'                                                          xi = .5,
#'                                                          w1 = 1,
#'                                                          w2 = 1,
#'                                                          w = 1,
#'                                                          lambda.beta = 1,
#'                                                          lambda.f = 1,
#'                                                          eta.beta = 1,
#'                                                          eta.f = 1,
#'                                                          n.lambda = 5,
#'                                                          n.eta = 5,
#'                                                          lambda.min.ratio = .01,
#'                                                          tol = 1e-3,
#'                                                          maxiter = 500,
#'                                                          report.prog = TRUE)
#'
#' plot_semipaddgt2pop_grid(semipadd2pop_logreg_grid.out,
#'                          true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#'                                                f2 = semipadd2pop_logreg_data$f2,
#'                                                X1 = semipadd2pop_logreg_data$X1,
#'                                                X2 = semipadd2pop_logreg_data$X2)
#' )
#' @export
semipadd2pop_grid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,response,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,eta.min.ratio = 0.001, eta.max.ratio = 10,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
  
  # get group lasso estimators over a grid of lambda and eta values
  if( response == "continuous"){
    
    grouplasso2pop_grid.out <- grouplasso2pop_linreg_grid(Y1 = Y1,
                                                          X1 = grouplasso2pop_inputs$DD1.tilde,
                                                          groups1 = grouplasso2pop_inputs$groups1,
                                                          Y2 = Y2,
                                                          X2 = grouplasso2pop_inputs$DD2.tilde,
                                                          groups2 = grouplasso2pop_inputs$groups2,
                                                          rho1 = rho1,
                                                          rho2 = rho2,
                                                          n.lambda = n.lambda,
                                                          n.eta = n.eta,
                                                          lambda.min.ratio = lambda.min.ratio,
                                                          lambda.max.ratio = lambda.max.ratio,
                                                          eta.min.ratio = eta.min.ratio,
                                                          eta.max.ratio = eta.max.ratio,
                                                          w1 = grouplasso2pop_inputs$w1,
                                                          w2 = grouplasso2pop_inputs$w2,
                                                          w = grouplasso2pop_inputs$w,
                                                          AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                          AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                          Com = grouplasso2pop_inputs$Com,
                                                          tol = tol,
                                                          maxiter = maxiter,
                                                          report.prog = report.prog)
                                                                 
  } else if( response == "binary"){
  
  grouplasso2pop_grid.out <- grouplasso2pop_logreg_grid(Y1 = Y1,
                                                        X1 = grouplasso2pop_inputs$DD1.tilde,
                                                        groups1 = grouplasso2pop_inputs$groups1,
                                                        Y2 = Y2,
                                                        X2 = grouplasso2pop_inputs$DD2.tilde,
                                                        groups2 = grouplasso2pop_inputs$groups2,
                                                        rho1 = rho1,
                                                        rho2 = rho2,
                                                        n.lambda = n.lambda,
                                                        n.eta = n.eta,
                                                        lambda.min.ratio = lambda.min.ratio,
                                                        lambda.max.ratio = lambda.max.ratio,
                                                        eta.min.ratio = eta.min.ratio,
                                                        eta.max.ratio = eta.max.ratio,
                                                        w1 = grouplasso2pop_inputs$w1,
                                                        w2 = grouplasso2pop_inputs$w2,
                                                        w = grouplasso2pop_inputs$w,
                                                        AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                        AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                        Com = grouplasso2pop_inputs$Com,
                                                        tol = tol,
                                                        maxiter = maxiter,
                                                        report.prog = report.prog)
  
  } else if(response == "gt"){
    
    grouplasso2pop_grid.out <- grouplasso2pop_gt_grid(Y1 = Y1$I,
                                                      Z1 = Y1$A,
                                                      Se1 = Y1$Se,
                                                      Sp1 = Y1$Sp,
                                                      E.approx1 = Y1$E.approx,
                                                      X1 = grouplasso2pop_inputs$DD1.tilde,
                                                      groups1 = grouplasso2pop_inputs$groups1,
                                                      Y2 = Y2$I,
                                                      Z2 = Y2$A,
                                                      Se2 = Y2$Se,
                                                      Sp2 = Y2$Sp,
                                                      E.approx2 = Y2$E.approx,
                                                      X2 = grouplasso2pop_inputs$DD2.tilde,
                                                      groups2 = grouplasso2pop_inputs$groups2,
                                                      rho1 = rho1,
                                                      rho2 = rho2,
                                                      n.lambda = n.lambda,
                                                      n.eta = n.eta,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      lambda.max.ratio = lambda.max.ratio,
                                                      eta.min.ratio = eta.min.ratio,
                                                      eta.max.ratio = eta.max.ratio,
                                                      w1 = grouplasso2pop_inputs$w1,
                                                      w2 = grouplasso2pop_inputs$w2,
                                                      w = grouplasso2pop_inputs$w,
                                                      AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                      AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                      Com = grouplasso2pop_inputs$Com,
                                                      tol = tol,
                                                      maxiter = maxiter,
                                                      report.prog = report.prog)
    
  } else {
    
    stop("invalid response type")
    
  }
  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

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
                                                              b1 = grouplasso2pop_grid.out$b1[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_grid.out$b2[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
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
                  lambda.seq = grouplasso2pop_grid.out$lambda.seq,
                  eta.seq = grouplasso2pop_grid.out$eta.seq,
                  iterations = grouplasso2pop_grid.out$iterations)

  class(output) <- "semipadd2pop_grid"

  return(output)

}

#' Compute semiparametric regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters
#'
#' @param Y1 the response data for data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the response data for data set 2
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
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
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param eta.min.ratio ratio of the smallest to largest value in the sequence of eta values
#' @param eta.max.ratio controls the largest value of eta in the eta sequence
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param report.prog a logical indicating whether the progress of the algorithm should be printed to the console
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' data <- get_semipadd2pop_data(n1 = 1000,n2 = 800,response = "gt")
#' 
#' semipadd2pop_cv.out <- semipadd2pop_cv(Y1 = data$Y1,
#'                                        X1 = data$X1,
#'                                        nonparm1 = data$nonparm1,
#'                                        Y2 = data$Y2,
#'                                        X2 = data$X2,
#'                                        nonparm2 = data$nonparm2,
#'                                        response = "gt",
#'                                        rho1 = 2,
#'                                        rho2 = 1,
#'                                        w1 = 1,
#'                                        w2 = 1,
#'                                        w = 1,
#'                                        nCom = data$nCom,
#'                                        d1 = data$nonparm1*25,
#'                                        d2 = data$nonparm2*15,
#'                                        xi = .5,
#'                                        n.lambda = 5,
#'                                        n.eta = 5,
#'                                        lambda.min.ratio = 0.01,
#'                                        lambda.max.ratio = 0.10,
#'                                        n.folds = 5,
#'                                        lambda.beta = 1,
#'                                        lambda.f = 1,
#'                                        eta.beta = 1,
#'                                        eta.f = 1,
#'                                        tol = 1e-3,
#'                                        maxiter = 1000,
#'                                        report.prog = FALSE)
#' 
#' plot_semipadd2pop_cv(semipadd2pop_cv.out) 
#' @export
semipadd2pop_cv <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,response,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,eta.min.ratio = 0.001, eta.max.ratio = 10,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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

  # get cv choice of group lasso tuning parameters
  if(response == "continuous"){
    
    grouplasso2pop_cv.out <- grouplasso2pop_linreg_cv(Y1 = Y1,
                                                      X1 = grouplasso2pop_inputs$DD1.tilde,
                                                      groups1 = grouplasso2pop_inputs$groups1,
                                                      Y2 = Y2,
                                                      X2 = grouplasso2pop_inputs$DD2.tilde,
                                                      groups2 = grouplasso2pop_inputs$groups2,
                                                      rho1 = rho1,
                                                      rho2 = rho2,
                                                      n.lambda = n.lambda,
                                                      n.eta = n.eta,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      lambda.max.ratio = lambda.max.ratio,
                                                      eta.min.ratio = eta.min.ratio,
                                                      eta.max.ratio = eta.max.ratio,
                                                      n.folds = n.folds,
                                                      w1 = grouplasso2pop_inputs$w1,
                                                      w2 = grouplasso2pop_inputs$w2,
                                                      w = grouplasso2pop_inputs$w,
                                                      AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                      AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                      Com = grouplasso2pop_inputs$Com,
                                                      tol = tol,
                                                      maxiter = maxiter,
                                                      report.prog = report.prog)
                                                             
  } else if(response == "binary"){
    
    
    grouplasso2pop_cv.out <- grouplasso2pop_logreg_cv(Y1 = Y1,
                                                      X1 = grouplasso2pop_inputs$DD1.tilde,
                                                      groups1 = grouplasso2pop_inputs$groups1,
                                                      Y2 = Y2,
                                                      X2 = grouplasso2pop_inputs$DD2.tilde,
                                                      groups2 = grouplasso2pop_inputs$groups2,
                                                      rho1 = rho1,
                                                      rho2 = rho2,
                                                      n.lambda = n.lambda,
                                                      n.eta = n.eta,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      lambda.max.ratio = lambda.max.ratio,
                                                      eta.min.ratio = eta.min.ratio,
                                                      eta.max.ratio = eta.max.ratio,
                                                      n.folds = n.folds,
                                                      w1 = grouplasso2pop_inputs$w1,
                                                      w2 = grouplasso2pop_inputs$w2,
                                                      w = grouplasso2pop_inputs$w,
                                                      AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                      AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                      Com = grouplasso2pop_inputs$Com,
                                                      tol = tol,
                                                      maxiter = maxiter,
                                                      report.prog = report.prog)
                                                             
    
  } else if(response == "gt"){
    
    grouplasso2pop_cv.out <- grouplasso2pop_gt_cv(Y1 = Y1$I,
                                                  Z1 = Y1$A,
                                                  Se1 = Y1$Se,
                                                  Sp1 = Y1$Sp,
                                                  E.approx1 = Y1$E.approx,
                                                  X1 = grouplasso2pop_inputs$DD1.tilde,
                                                  groups1 = grouplasso2pop_inputs$groups1,
                                                  Y2 = Y2$I,
                                                  Z2 = Y2$A,
                                                  Se2 = Y2$Se,
                                                  Sp2 = Y2$Sp,
                                                  E.approx2 = Y2$E.approx,
                                                  X2 = grouplasso2pop_inputs$DD2.tilde,
                                                  groups2 = grouplasso2pop_inputs$groups2,
                                                  rho1 = rho1,
                                                  rho2 = rho2,
                                                  n.lambda = n.lambda,
                                                  n.eta = n.eta,
                                                  lambda.min.ratio = lambda.min.ratio,
                                                  lambda.max.ratio = lambda.max.ratio,
                                                  eta.min.ratio = eta.min.ratio,
                                                  eta.max.ratio = eta.max.ratio,
                                                  n.folds = n.folds,
                                                  w1 = grouplasso2pop_inputs$w1,
                                                  w2 = grouplasso2pop_inputs$w2,
                                                  w = grouplasso2pop_inputs$w,
                                                  AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                  AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                  Com = grouplasso2pop_inputs$Com,
                                                  tol = tol,
                                                  maxiter = maxiter,
                                                  report.prog = report.prog)
  } else{
    
    stop("invalid response type")
    
  }
  
  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplasso2pop_cv.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_cv.out$b2.arr[,l,k])

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
                  beta1.hat = beta1.hat,
                  beta2.hat = beta2.hat,
                  nonparm1 = nonparm1,
                  nonparm2 = nonparm2,
                  rho1 = rho1,
                  rho2 = rho2,
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
                  lambda.seq = grouplasso2pop_cv.out$lambda.seq,
                  eta.seq = grouplasso2pop_cv.out$eta.seq,
                  which.lambda.cv = grouplasso2pop_cv.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_cv.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_cv.out$which.lambda.cv.under.zero.eta,
                  iterations = grouplasso2pop_cv.out$iterations)

  class(output) <- "sempadd2pop_cv"

  return(output)

}
 
#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#'
#' @param Y1 the response data for data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the response data for data set 2
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
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
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param eta.min.ratio ratio of the smallest to largest value in the sequence of eta values
#' @param eta.max.ratio controls the largest value of eta in the eta sequence
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param report.prog a logical indicating whether the progress of the algorithm should be printed to the console
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' data <- get_semipadd2pop_data(n1 = 1000,n2 = 800,response = "gt")
#' 
#' semipadd2pop_cv_adapt.out <- semipadd2pop_cv_adapt(Y1 = data$Y1,
#'                                                    X1 = data$X1,
#'                                                    nonparm1 = data$nonparm1,
#'                                                    Y2 = data$Y2,
#'                                                    X2 = data$X2,
#'                                                    nonparm2 = data$nonparm2,
#'                                                    response = "gt",
#'                                                    rho1 = 2,
#'                                                    rho2 = 1,
#'                                                    w1 = 1,
#'                                                    w2 = 1,
#'                                                    w = 1,
#'                                                    nCom = data$nCom,
#'                                                    d1 = 25,
#'                                                    d2 = 15,
#'                                                    xi = .5,
#'                                                    n.lambda = 5,
#'                                                    n.eta = 5,
#'                                                    lambda.min.ratio = 0.01,
#'                                                    lambda.max.ratio = 0.10,
#'                                                    n.folds = 5,
#'                                                    lambda.beta = 1,
#'                                                    lambda.f = 1,
#'                                                    eta.beta = 1,
#'                                                    eta.f = 1,
#'                                                    tol = 1e-3,
#'                                                    maxiter = 1000,
#'                                                    report.prog = TRUE)
#' 
#' plot_semipadd2pop_cv(semipadd2pop_cv_adapt.out) 
#' @export
semipadd2pop_cv_adapt <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,response,rho1,rho2,w1=1,w2=1,w=1,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio = .01,lambda.max.ratio=1,eta.min.ratio = 0.001, eta.max.ratio =10,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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

  # get group lasso estimators over a grid of lambda and eta values
  if(response == "continuous"){
    
    grouplasso2pop_cv_adapt.out <- grouplasso2pop_linreg_cv_adapt(Y1 = Y1,
                                                                  X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                  groups1 = grouplasso2pop_inputs$groups1,
                                                                  Y2 = Y2,
                                                                  X2 = grouplasso2pop_inputs$DD2.tilde,
                                                                  groups2 = grouplasso2pop_inputs$groups2,
                                                                  rho1 = rho1,
                                                                  rho2 = rho2,
                                                                  n.lambda = n.lambda,
                                                                  n.eta = n.eta,
                                                                  lambda.min.ratio = lambda.min.ratio,
                                                                  lambda.max.ratio = lambda.max.ratio,
                                                                  eta.min.ratio = eta.min.ratio,
                                                                  eta.max.ratio = eta.max.ratio,
                                                                  n.folds = n.folds,
                                                                  w1 = grouplasso2pop_inputs$w1,
                                                                  w2 = grouplasso2pop_inputs$w2,
                                                                  w = grouplasso2pop_inputs$w,
                                                                  AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                                  AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                                  Com = grouplasso2pop_inputs$Com,
                                                                  tol = tol,
                                                                  maxiter = maxiter,
                                                                  report.prog = report.prog) 
    
    
  } else if(response == "binary"){
    
    grouplasso2pop_cv_adapt.out <- grouplasso2pop_logreg_cv_adapt(Y1 = Y1,
                                                                  X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                  groups1 = grouplasso2pop_inputs$groups1,
                                                                  Y2 = Y2,
                                                                  X2 = grouplasso2pop_inputs$DD2.tilde,
                                                                  groups2 = grouplasso2pop_inputs$groups2,
                                                                  rho1 = rho1,
                                                                  rho2 = rho2,
                                                                  n.lambda = n.lambda,
                                                                  n.eta = n.eta,
                                                                  lambda.min.ratio = lambda.min.ratio,
                                                                  lambda.max.ratio = lambda.max.ratio,
                                                                  eta.min.ratio = eta.min.ratio,
                                                                  eta.max.ratio = eta.max.ratio,
                                                                  n.folds = n.folds,
                                                                  w1 = grouplasso2pop_inputs$w1,
                                                                  w2 = grouplasso2pop_inputs$w2,
                                                                  w = grouplasso2pop_inputs$w,
                                                                  AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                                  AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                                  Com = grouplasso2pop_inputs$Com,
                                                                  tol = tol,
                                                                  maxiter = maxiter,
                                                                  report.prog = report.prog) 
                                                                         
    
  } else if(response == "gt"){
    
    grouplasso2pop_cv_adapt.out <- grouplasso2pop_gt_cv_adapt(Y1 = Y1$I,
                                                              Z1 = Y1$A,
                                                              Se1 = Y1$Se,
                                                              Sp1 = Y1$Sp,
                                                              E.approx1 = Y1$E.approx,
                                                              X1 = grouplasso2pop_inputs$DD1.tilde,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              Y2 = Y2$I,
                                                              Z2 = Y2$A,
                                                              Se2 = Y2$Se,
                                                              Sp2 = Y2$Sp,
                                                              E.approx2 = Y2$E.approx,
                                                              X2 = grouplasso2pop_inputs$DD2.tilde,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              rho1 = rho1,
                                                              rho2 = rho2,
                                                              n.lambda = n.lambda,
                                                              n.eta = n.eta,
                                                              lambda.min.ratio = lambda.min.ratio,
                                                              lambda.max.ratio = lambda.max.ratio,
                                                              eta.min.ratio = eta.min.ratio,
                                                              eta.max.ratio = eta.max.ratio,
                                                              n.folds = n.folds,
                                                              w1 = grouplasso2pop_inputs$w1,
                                                              w2 = grouplasso2pop_inputs$w2,
                                                              w = grouplasso2pop_inputs$w,
                                                              AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                              AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                              Com = grouplasso2pop_inputs$Com,
                                                              tol = tol,
                                                              maxiter = maxiter,
                                                              report.prog = report.prog)
                                                  
    
  } else{
    
    stop("invalid response type")
    
  }

  # get matrices of the fitted functions evaluated at the design points
  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplasso2pop_cv_adapt.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_cv_adapt.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
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
                  lambda.seq = grouplasso2pop_cv_adapt.out$lambda.seq,
                  eta.seq = grouplasso2pop_cv_adapt.out$eta.seq,
                  lambda.initial.fit = grouplasso2pop_cv_adapt.out$lambda.initial.fit,
                  which.lambda.cv = grouplasso2pop_cv_adapt.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_cv_adapt.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_cv_adapt.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplasso2pop_cv_adapt.out$w1,
                  w2 = grouplasso2pop_cv_adapt.out$w2,
                  w = grouplasso2pop_cv_adapt.out$w,
                  iterations = grouplasso2pop_cv_adapt.out$iterations)

  class(output) <- "semipadd2pop_cv"

  return(output)

}

#' Compute semiparametric regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters from a user-specified grid after an adaptive step
#'
#' @param Y1 the response data for data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the response data for data set 2
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{XX1} and \code{XX2} after the column of ones corresponding to the intercept.
#' @param d the dimension of the B-spline basis to be used when fitting the nonparametric effects
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.seq the sequence of lambda values
#' @param eta.seq the sequence of eta values
#' @param lambda.initial.fit the lambda value to use for the initial estimator used for making the weights for the adaptive step
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501, n2 = 604)
#'
#' semipadd2pop_logreg_cv_adapt.out <- semipadd2pop_logreg_cv_adapt( Y1 = semipadd2pop_logreg_data$Y1,
#'                                                                   X1 = semipadd2pop_logreg_data$X1,
#'                                                                   nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                                   Y2 = semipadd2pop_logreg_data$Y2,
#'                                                                   X2 = semipadd2pop_logreg_data$X2,
#'                                                                   nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                                   w1 = 1,
#'                                                                   w2 = 1,
#'                                                                   w = 1,
#'                                                                   nCom = semipadd2pop_logreg_data$nCom,
#'                                                                   d1 = semipadd2pop_logreg_data$nonparm1*25,
#'                                                                   d2 = semipadd2pop_logreg_data$nonparm2*15,
#'                                                                   xi = .5,
#'                                                                   n.lambda = 5,
#'                                                                   n.eta = 5,
#'                                                                   lambda.min.ratio = .01,
#'                                                                   n.folds = 5,
#'                                                                   lambda.beta = 1,
#'                                                                   lambda.f = 1,
#'                                                                   eta.beta = 1,
#'                                                                   eta.f = 1,
#'                                                                   tol = 1e-3,
#'                                                                   maxiter = 1000,
#'                                                                   report.prog = FALSE)
#'
#' semipadd2pop_logreg_cv_adapt_fixedgrid.out <- semipadd2pop_logreg_cv_adapt_fixedgrid(Y1 = semipadd2pop_logreg_data$Y1,
#'                                                                                      X1 = semipadd2pop_logreg_data$X1,
#'                                                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                                                      Y2 = semipadd2pop_logreg_data$Y2,
#'                                                                                      X2 = semipadd2pop_logreg_data$X2,
#'                                                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                                                      rho1 = 2,
#'                                                                                      rho2 = 1,
#'                                                                                      nCom = semipadd2pop_logreg_data$nCom,
#'                                                                                      d = 20,
#'                                                                                      xi = .5,
#'                                                                                      w1 = 1,
#'                                                                                      w2 = 1,
#'                                                                                      w = 1,
#'                                                                                      lambda.seq = semipadd2pop_logreg_cv_adapt.out$lambda.seq,
#'                                                                                      eta.seq = semipadd2pop_logreg_cv_adapt.out$eta.seq,
#'                                                                                      n.folds = 5,
#'                                                                                      lambda.beta = 1,
#'                                                                                      lambda.f = 1,
#'                                                                                      eta.beta = 1,
#'                                                                                      eta.f = 1,
#'                                                                                      tol = 1e-2,
#'                                                                                      maxiter = 500)
#'
#' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv_adapt_fixedgrid.out,
#'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#'                                              f2 = semipadd2pop_logreg_data$f2,
#'                                              X1 = semipadd2pop_logreg_data$X1,
#'                                              X2 = semipadd2pop_logreg_data$X2)
#' )
#' @export
semipadd2pop_cv_adapt_fixedgrid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,response,rho1,rho2,nCom,d1,d2,xi,w1,w2,w,lambda.seq,eta.seq,lambda.initial.fit,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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

  # get group lasso estimators over a grid of lambda and eta values
  if(response == "continuous"){
    
    grouplasso2pop_cv_adapt_fixedgrid.out <- grouplasso2pop_linreg_cv_adapt_fixedgrid(Y1 = Y1,
                                                                                      X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                                      groups1 = grouplasso2pop_inputs$groups1,
                                                                                      Y2 = Y2,
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
                                                                                      tol = tol,
                                                                                      maxiter = maxiter,
                                                                                      report.prog = report.prog)
    
  } else if(response == "binary"){
    
    
    grouplasso2pop_cv_adapt_fixedgrid.out <- grouplasso2pop_logreg_cv_adapt_fixedgrid(Y1 = Y1,
                                                                                      X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                                      groups1 = grouplasso2pop_inputs$groups1,
                                                                                      Y2 = Y2,
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
                                                                                      tol = tol,
                                                                                      maxiter = maxiter,
                                                                                      report.prog = report.prog)
                                                                                             
  } else if(response == "gt"){
    
    grouplasso2pop_cv_adapt_fixedgrid.out <- grouplasso2pop_gt_cv_adapt_fixedgrid(Y1 = Y1$I,
                                                                                  Z1 = Y1$A,
                                                                                  Se1 = Y1$Se,
                                                                                  Sp1 = Y1$Sp,
                                                                                  E.approx1 = Y1$E.approx,
                                                                                  X1 = grouplasso2pop_inputs$DD1.tilde,
                                                                                  groups1 = grouplasso2pop_inputs$groups1,
                                                                                  Y2 = Y2$I,
                                                                                  Z2 = Y2$A,
                                                                                  Se2 = Y2$Se,
                                                                                  Sp2 = Y2$Sp,
                                                                                  E.approx2 = Y2$E.approx,
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
                                                                                  tol = tol,
                                                                                  maxiter = maxiter,
                                                                                  report.prog = report.prog)
                                                                                      
  } else{
    
    
    stop("invalid response type")
    
  }

  which.lambda.cv <- grouplasso2pop_cv_adapt_fixedgrid.out$which.lambda.cv
  which.eta.cv <- grouplasso2pop_cv_adapt_fixedgrid.out$which.eta.cv

  # get matrices of the fitted functions evaluated at the design points

  n.lambda <- length(lambda.seq)
  n.eta <- length(eta.seq)

  f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
  f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))

  beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))

  f1.hat <- vector("list",n.lambda)
  f2.hat <- vector("list",n.lambda)

  for(l in 1:n.lambda)
    for(k in 1:n.eta)
    {

      semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                              nonparm1 = nonparm1,
                                                              groups1 = grouplasso2pop_inputs$groups1,
                                                              knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                              emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                              QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                              b1 = grouplasso2pop_cv_adapt_fixedgrid.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_cv_adapt_fixedgrid.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

    }

  # prepare output
  output <- list( f1.hat = f1.hat,
                  f2.hat = f2.hat,
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
                  which.lambda.cv = grouplasso2pop_cv_adapt_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_cv_adapt_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplasso2pop_cv_adapt_fixedgrid.out$w1,
                  w2 = grouplasso2pop_cv_adapt_fixedgrid.out$w2,
                  w = grouplasso2pop_cv_adapt_fixedgrid.out$w,
                  iterations = grouplasso2pop_cv_adapt_fixedgrid.out$iterations)

  class(output) <- "semipadd2pop_gt_cv"

  return(output)

}