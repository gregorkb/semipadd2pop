
#' #' Prepare inputs for grouplasso2pop function when using it to fit a semiparametric model
#' #'
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @export
#' semipadd2pop_to_grouplasso2pop <- function(X1,nonparm1,X2,nonparm2,nCom,d1,d2,ComInt=NULL,Int1=NULL,Int2=NULL,xi,w1=1,w2=1,w=1,w1_int=1,w2_int=1,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1)
#' {
#'   
#'   ww1 <- ((1-nonparm1) + nonparm1 * lambda.f/lambda.beta) * w1
#'   ww2 <- ((1-nonparm2) + nonparm2 * lambda.f/lambda.beta) * w2
#'   ww1[1] <- 0 # do not penalize intercept
#'   ww2[1] <- 0 # do not penalize intercept
#'   ww <- ((1-nonparm1) + nonparm1 * eta.f/eta.beta) * w
#'   Com <- 1 + 1:nCom
#'   ww <- ww[1:max(Com)]
#'   
#'   n1 <- nrow(X1)
#'   n2 <- nrow(X2)
#'   
#'   pp1 <- ncol(X1)
#'   pp2 <- ncol(X2)
#'   
#'   ##### For first data set
#'   DD1.tilde <- matrix(NA,n1,0)
#'   groups1 <- numeric()
#'   AA1.tilde <- vector("list",length=pp1)
#'   QQ1.inv <- vector("list",length=pp1)
#'   knots.list1 <- vector("list",length=pp1)
#'   emp.cent1 <- vector("list",length=pp1)
#'   
#'   if( length(d1) == 1 ){
#'     
#'     d1 <- rep(d1,pp1)
#'     
#'   }
#'   
#'   if( length(d2) == 1 ){
#'     
#'     d2 <- rep(d2,pp1)
#'     
#'   }
#'   
#'   for( j in 1:pp1 )
#'   {
#'     
#'     if(nonparm1[j] == 0){
#'       
#'       DD1.tilde <- cbind(DD1.tilde,X1[,j])
#'       AA1.tilde[[j]] <- as.matrix(1)
#'       groups1 <- c(groups1,j)
#'       
#'     } else {
#'       
#'       if(d1[j] < 0){
#'         
#'         int.knots <- seq(min(X1[,j]), max(X1[,j]), length = -d1[j] - 2 + 1 )
#'         d1[j] <- -d1[j]
#'         
#'       } else {
#'         
#'         int.knots <- quantile(X1[,j],seq(0,1,length=d1[j]-2+1)) # add one, so that one can be removed after centering to restore full-rank.
#'         
#'       }
#'       
#'       boundary.knots <- range(int.knots)
#'       all.knots <- sort(c(rep(boundary.knots,3),int.knots))
#'       knots.list1[[j]] <- all.knots
#'       
#'       B1j <- spline.des(all.knots,X1[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
#'       emp.cent1[[j]] <- apply(B1j,2,mean)
#'       B1j.cent <- B1j - matrix(emp.cent1[[j]],n1,d1[j],byrow=TRUE)
#'       
#'       # construct matrix in which l2 norm of function is a quadratic form
#'       M <- t(B1j.cent) %*% B1j.cent / n1
#'       
#'       # construct matrix in which 2nd derivative penalty is a quadratic form
#'       R <- matrix(NA,d1[j]+1,d1[j]+1)
#'       dsq_bspline.mat <- spline.des(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)$design
#'       for(k in 1:(d1[j]+1))
#'         for(l in 1:(d1[j]+1))
#'         {
#'           
#'           pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
#'           h <- diff(int.knots)
#'           R[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.
#'           
#'         }
#'       
#'       Q <-  chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
#'       Q.inv <- solve(Q)
#'       
#'       QQ1.inv[[j]] <- Q.inv
#'       
#'       # construct transformed basis function matrices
#'       DD1.tilde <- cbind(DD1.tilde, B1j.cent %*% Q.inv)
#'       
#'       if(j %in% Com){  # Make the AA matrices based on a common centering between the data sets... since we are not really interested in the centering but in the shape.
#'         # We penalize the shapes to be similar, not the intercepts!!!!!!  I don't know if this will address the problem I am seeing though....
#'         
#'         X1j.srt <- sort(X1[,j])
#'         X2j.srt <- sort(X2[,j])
#'         
#'         Xj.union <- sort(c(X1[,j],X2[,j]))
#'         int.ind <- which( max(X1j.srt[1],X2j.srt[1]) < Xj.union & Xj.union < min(X1j.srt[n1],X2j.srt[n2]) )
#'         
#'         Xj.int <- Xj.union[int.ind]
#'         
#'         AA1j <- spline.des(all.knots,Xj.int,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove same one
#'         AA1j.cent <- AA1j - matrix(apply(AA1j,2,mean),length(int.ind),d1[j],byrow=TRUE)
#'         AA1.tilde[[j]] <- AA1j.cent %*% Q.inv
#'         
#'       }
#'       
#'       groups1 <- c(groups1,rep(j,d1[j]))
#'       
#'     }
#'     
#'   }
#'   
#'   ##### For second data set
#'   DD2.tilde <- matrix(NA,n2,0)
#'   groups2 <- numeric()
#'   AA2.tilde <- vector("list",length=pp2)
#'   QQ2.inv <- vector("list",length=pp2)
#'   knots.list2 <- vector("list",length=pp2)
#'   emp.cent2 <- vector("list",length=pp2)
#'   
#'   for( j in 1:pp2 )
#'   {
#'     
#'     if(nonparm2[j] == 0){
#'       
#'       DD2.tilde <- cbind(DD2.tilde,X2[,j])
#'       AA2.tilde[[j]] <- as.matrix(1)
#'       groups2 <- c(groups2,j)
#'       
#'     } else {
#'       
#'       if(d2[j] < 0){
#'         
#'         int.knots <- seq(min(X2[,j]), max(X2[,j]), length = -d2[j] - 2 + 1 )
#'         d2[j] <- -d2[j]
#'         
#'       } else {
#'         
#'         int.knots <- quantile(X2[,j],seq(0,1,length=d2[j]-2+1)) # add one, so that one can be removed after centering to restore full-rank.
#'         
#'       }
#'       
#'       boundary.knots <- range(int.knots)
#'       all.knots <- sort(c(rep(boundary.knots,3),int.knots))
#'       knots.list2[[j]] <- all.knots
#'       
#'       B2j <- spline.des(all.knots,X2[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
#'       emp.cent2[[j]] <- apply(B2j,2,mean)
#'       B2j.cent <- B2j - matrix(emp.cent2[[j]],n2,d2[j],byrow=TRUE)
#'       
#'       # construct matrix in which l2 norm of function is a quadratic form
#'       M <- t(B2j.cent) %*% B2j.cent / n2
#'       
#'       # construct matrix in which 2nd derivative penalty is a quadratic form
#'       R <- matrix(NA,d2[j]+1,d2[j]+1)
#'       dsq_bspline.mat <- spline.des(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)$design
#'       for(k in 1:(d2[j]+1))
#'         for(l in 1:(d2[j]+1))
#'         {
#'           
#'           pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
#'           h <- diff(int.knots)
#'           R[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.
#'           
#'         }
#'       
#'       Q <-  chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
#'       Q.inv <- solve(Q)
#'       
#'       QQ2.inv[[j]] <- Q.inv
#'       
#'       # construct transformed basis function matrices
#'       DD2.tilde <- cbind(DD2.tilde, B2j.cent %*% Q.inv)
#'       
#'       if(j %in% Com){  # Make the AA matrices based on a common centering between the data sets... since we are not really interested in the centering but in the shape.
#'         # We penalize the shapes to be similar, not the intercepts!!!!!!  I don't know if this will address the problem I am seeing though....
#'         
#'         X1j.srt <- sort(X1[,j])
#'         X2j.srt <- sort(X2[,j])
#'         
#'         Xj.union <- sort(c(X1[,j],X2[,j]))
#'         int.ind <- which( max(X1j.srt[1],X2j.srt[1]) < Xj.union & Xj.union < min(X1j.srt[n1],X2j.srt[n2]) )
#'         
#'         Xj.int <- Xj.union[int.ind]
#'         
#'         AA2j <- spline.des(all.knots,Xj.int,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove same one
#'         AA2j.cent <- AA2j - matrix(apply(AA2j,2,mean),length(int.ind),d2[j],byrow=TRUE)
#'         AA2.tilde[[j]] <- AA2j.cent %*% Q.inv
#'         
#'       }
#'       
#'       groups2 <- c(groups2,rep(j,d2[j]))
#'       
#'     }
#'     
#'   }
#'   
#'   output <- list( DD1.tilde = DD1.tilde,
#'                   DD2.tilde = DD2.tilde,
#'                   groups1 = groups1,
#'                   groups2 = groups2,
#'                   knots.list1 = knots.list1,
#'                   knots.list2 = knots.list2,
#'                   QQ1.inv = QQ1.inv,
#'                   QQ2.inv = QQ2.inv,
#'                   emp.cent1 = emp.cent1,
#'                   emp.cent2 = emp.cent2,
#'                   AA1.tilde = AA1.tilde,
#'                   AA2.tilde = AA2.tilde,
#'                   lambda = lambda.beta,
#'                   eta = eta.beta,
#'                   w1 = ww1,
#'                   w2 = ww2,
#'                   w = ww,
#'                   Com = Com)
#'   
#'   return(output)
#'   
#' }
#' 
#' #' Convert output from grouplasso2pop to the fitted functions of the semi-parametric additive model
#' #'
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param groups1 a vector indicating to which group the entries of the coefficient vector \code{b1} belong
#' #' @param knots.list1 a list of vectors with the knot locations for nonparametric effects for data set 1
#' #' @param emp.cent1 a list of vectors of the empirical basis function centerings for data set 1
#' #' @param QQ1.inv the matrix with which to back-transform the group lasso coefficients for data set 1
#' #' @param b1 the group lasso coefficients for data set 1
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param groups2 a vector indicating to which group the entries of the coefficient vector \code{b2} belong
#' #' @param knots.list2 a list of vectors with the knot locations for nonparametric effects for data set 2
#' #' @param emp.cent2 a list of vectors of the empirical basis function centerings for data set 2
#' #' @param QQ2.inv the matrix with which to back-transform the group lasso coefficients for data set 2
#' #' @param b2 the group lasso coefficients for data set 2
#' #' @export
#' grouplasso2pop_to_semipadd2pop <- function(X1,nonparm1,groups1,knots.list1,emp.cent1,QQ1.inv,b1,X2,nonparm2,groups2,knots.list2,emp.cent2,QQ2.inv,b2)
#' {
#'   
#'   n1 <- nrow(X1)
#'   n2 <- nrow(X2)
#'   
#'   pp1 <- ncol(X1)
#'   pp2 <- ncol(X2)
#'   
#'   # store fitted functions on data set 1 in a list
#'   f1.hat <- vector("list",pp1)
#'   f1.hat.design <- matrix(0,n1,pp1)
#'   
#'   f1.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b1[1])," }")))
#'   f1.hat.design[,1] <- b1[1]
#'   
#'   
#'   beta1.hat <- rep(NA,pp1)
#'   beta1.hat[1] <- b1[1]
#'   for(j in 2:pp1)
#'   {
#'     
#'     if(nonparm1[j] == 0)
#'     {
#'       
#'       ind <- which(groups1 == j)
#'       f1.hat[[j]] <- eval( parse( text = paste("function(x){ x * ",paste(b1[ind])," }")))
#'       beta1.hat[j] <- b1[ind]
#'       
#'     } else {
#'       
#'       ind <- which(groups1 == j)
#'       d <- length(ind)
#'       
#'       Gamma1.hat <- QQ1.inv[[j]] %*% b1[ind]
#'       f1.hat[[j]] <- eval(parse(text = paste("function(x)
#'                                              {
#' 
#'                                              x <- round(x,10)
#'                                              x.mat <- spline.des(",paste("c(",paste(round(knots.list1[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
#'                                              x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent1[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
#'                                              f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma1.hat,collapse=","),")",sep=""),")
#' 
#'                                              return(f.hat)
#' 
#'                                              }"
#'       )))
#'       
#'     }
#'     
#'     f1.hat.design[,j] <- f1.hat[[j]](X1[,j])
#'     
#'   }
#'   
#'   # store fitted functions on data set 2 in a list
#'   f2.hat <- vector("list",pp2)
#'   f2.hat.design <- matrix(0,n2,pp2)
#'   
#'   f2.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b2[1])," }")))
#'   f2.hat.design[,1] <- b2[1]
#'   
#'   beta2.hat <- rep(NA,pp2)
#'   beta2.hat[1] <- b2[1]
#'   for(j in 2:pp2)
#'   {
#'     
#'     if(nonparm2[j] == 0)
#'     {
#'       
#'       ind <- which(groups2 == j)
#'       f2.hat[[j]] <- eval( parse( text = paste("function(x){ x * ",paste(b2[ind])," }")))
#'       beta2.hat[j] <- b2[ind]
#'       
#'     } else {
#'       
#'       ind <- which(groups2 == j)
#'       d <- length(ind)
#'       
#'       Gamma2.hat <- QQ2.inv[[j]] %*% b2[ind]
#'       f2.hat[[j]] <- eval(parse(text=paste("function(x)
#'                                            {
#' 
#'                                            x <- round(x,10)
#'                                            x.mat <- spline.des(",paste("c(",paste(round(knots.list2[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
#'                                            x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent2[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
#'                                            f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma2.hat,collapse=","),")",sep=""),")
#' 
#'                                            return(f.hat)
#' 
#'                                            }"
#'       )))
#'       
#'     }
#'     
#'     f2.hat.design[,j] <- f2.hat[[j]](X2[,j])
#'     
#'   }
#'   
#'   output <- list(f1.hat = f1.hat,
#'                  f1.hat.design = f1.hat.design,
#'                  f2.hat = f2.hat,
#'                  f2.hat.design = f2.hat.design,
#'                  beta1.hat = beta1.hat,
#'                  beta2.hat = beta2.hat)
#'   
#' }
#' 
#' 
#' 
#' #' Compute semiparametric continuous-response regression model with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 the response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
#' #'
#' #' @examples
#' #' semipadd2pop_linreg_data <- get_semipadd2pop_linreg_data(n1 = 501, n2 = 604)
#' #'
#' #' semipadd2pop_linreg.out <- semipadd2pop_linreg(Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                X1 = semipadd2pop_linreg_data$X1,
#' #'                                                nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                X2 = semipadd2pop_linreg_data$X2,
#' #'                                                nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                w1=1,
#' #'                                                w2=1,
#' #'                                                w=1,
#' #'                                                nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                xi=.5,
#' #'                                                lambda.beta=1,
#' #'                                                lambda.f=1,
#' #'                                                eta.beta=1,
#' #'                                                eta.f=1,
#' #'                                                tol=1e-3,
#' #'                                                maxiter=500,
#' #'                                                plot_obj=FALSE)
#' #'
#' #' plot_semipaddgt2pop(semipadd2pop_linreg.out,
#' #'                     true.functions=list(f1 = semipadd2pop_linreg_data$f1,
#' #'                                         f2 = semipadd2pop_linreg_data$f2,
#' #'                                         X1 = semipadd2pop_linreg_data$X1,
#' #'                                         X2 = semipadd2pop_linreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_linreg <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,w1=1,w2=1,w=1,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,tol=1e-4,maxiter=500,plot_obj=FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators
#'   grouplasso2pop.out <- grouplasso2pop_linreg( Y1 = Y1,
#'                                                X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                Y2 = Y2,
#'                                                X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                lambda = grouplasso2pop_inputs$lambda,
#'                                                eta = grouplasso2pop_inputs$eta,
#'                                                w1 = grouplasso2pop_inputs$w1,
#'                                                w2 = grouplasso2pop_inputs$w2,
#'                                                w = grouplasso2pop_inputs$w,
#'                                                AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                Com = grouplasso2pop_inputs$Com,
#'                                                tol=tol,
#'                                                maxiter=maxiter,
#'                                                plot_obj=plot_obj)
#'   
#'   # construct fitted functions from grouplasso2pop output
#'   semipadd2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                         nonparm1 = nonparm1,
#'                                                         groups1 = grouplasso2pop_inputs$groups1,
#'                                                         knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                         emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                         QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                         b1 = grouplasso2pop.out$beta1.hat,
#'                                                         X2 = X2,
#'                                                         nonparm2 = nonparm2,
#'                                                         groups2 = grouplasso2pop_inputs$groups2,
#'                                                         knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                         emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                         QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                         b2 = grouplasso2pop.out$beta2.hat)
#'   
#'   # collect output
#'   output <- list( f1.hat = semipadd2pop_fitted$f1.hat,
#'                   f2.hat = semipadd2pop_fitted$f2.hat,
#'                   f1.hat.design = semipadd2pop_fitted$f1.hat.design,
#'                   f2.hat.design = semipadd2pop_fitted$f2.hat.design,
#'                   beta1.hat = semipadd2pop_fitted$beta1.hat,
#'                   beta2.hat = semipadd2pop_fitted$beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com)
#'   
#'   return(output)
#'   
#' }
#' 
#' #' Compute semiparametric continuous-response regression model with 2 data sets while penalizing dissimilarity over a grid of tuning parameter values
#' #'
#' #' @param Y1 the continuous response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the continuous response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_linreg_data <- get_semipadd2pop_linreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_linreg_grid.out <- semipadd2pop_linreg_grid(Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                          X1 = semipadd2pop_linreg_data$X1,
#' #'                                                          nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                          Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                          X2 = semipadd2pop_linreg_data$X2,
#' #'                                                          nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                          nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                          d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                          d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                          xi = .5,
#' #'                                                          w1 = 1,
#' #'                                                          w2 = 1,
#' #'                                                          w = 1,
#' #'                                                          lambda.beta = 1,
#' #'                                                          lambda.f = 1,
#' #'                                                          eta.beta = 1,
#' #'                                                          eta.f = 1,
#' #'                                                          n.lambda = 5,
#' #'                                                          n.eta = 5,
#' #'                                                          lambda.min.ratio = .01,
#' #'                                                          tol = 1e-3,
#' #'                                                          maxiter = 500,
#' #'                                                          report.prog = TRUE)
#' #'
#' #' plot_semipaddgt2pop_grid(semipadd2pop_linreg_grid.out,
#' #'                          true.functions = list(f1 = semipadd2pop_linreg_data$f1,
#' #'                                                f2 = semipadd2pop_linreg_data$f2,
#' #'                                                X1 = semipadd2pop_linreg_data$X1,
#' #'                                                X2 = semipadd2pop_linreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_linreg_grid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_linreg_grid.out <- grouplasso2pop_linreg_grid(Y1 = Y1,
#'                                                                X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                Y2 = Y2,
#'                                                                X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                n.lambda,
#'                                                                n.eta,
#'                                                                lambda.min.ratio,
#'                                                                w1 = grouplasso2pop_inputs$w1,
#'                                                                w2 = grouplasso2pop_inputs$w2,
#'                                                                w = grouplasso2pop_inputs$w,
#'                                                                AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                Com = grouplasso2pop_inputs$Com,
#'                                                                tol=tol,
#'                                                                maxiter=maxiter,
#'                                                                report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1,
#'                                                               nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_linreg_grid.out$b1[,l,k],
#'                                                               X2,
#'                                                               nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_linreg_grid.out$b2[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   lambda.seq = grouplasso2pop_linreg_grid.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_linreg_grid.out$eta.seq,
#'                   iterations = grouplasso2pop_linreg_grid.out$iterations)
#'   
#'   
#'   class(output) <- "semipaddgt2pop_grid"
#'   
#'   return(output)
#'   
#' }
#' 
#' #' Compute semiparametric continuous-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters
#' #'
#' #' @param Y1 the continuous response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the continuous response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_linreg_data <- get_semipadd2pop_linreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_linreg_cv.out <- semipadd2pop_linreg_cv(Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                      X1 = semipadd2pop_linreg_data$X1,
#' #'                                                      nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                      Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                      X2 = semipadd2pop_linreg_data$X2,
#' #'                                                      nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                      w1 = 1,
#' #'                                                      w2 = 1,
#' #'                                                      w = 1,
#' #'                                                      nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                      d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                      d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                      xi = .5,
#' #'                                                      n.lambda = 5,
#' #'                                                      n.eta = 5,
#' #'                                                      lambda.min.ratio = .01,
#' #'                                                      n.folds = 5,
#' #'                                                      lambda.beta = 1,
#' #'                                                      lambda.f = 1,
#' #'                                                      eta.beta = 1,
#' #'                                                      eta.f = 1,
#' #'                                                      tol = 1e-3,
#' #'                                                      maxiter = 1000,
#' #'                                                      report.prog = FALSE)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_linreg_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_linreg_data$f1,
#' #'                                              f2 = semipadd2pop_linreg_data$f2,
#' #'                                              X1 = semipadd2pop_linreg_data$X1,
#' #'                                              X2 = semipadd2pop_linreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_linreg_cv <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_linreg_cv.out <- grouplasso2pop_linreg_cv(Y1 = Y1,
#'                                                            X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                            groups1 = grouplasso2pop_inputs$groups1,
#'                                                            Y2 = Y2,
#'                                                            X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                            groups2 = grouplasso2pop_inputs$groups2,
#'                                                            n.lambda = n.lambda,
#'                                                            n.eta = n.eta,
#'                                                            lambda.min.ratio = lambda.min.ratio,
#'                                                            n.folds = n.folds,
#'                                                            w1 = grouplasso2pop_inputs$w1,
#'                                                            w2 = grouplasso2pop_inputs$w2,
#'                                                            w = grouplasso2pop_inputs$w,
#'                                                            AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                            AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                            Com = grouplasso2pop_inputs$Com,
#'                                                            tol=tol,
#'                                                            maxiter=maxiter,
#'                                                            report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1,
#'                                                               nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_linreg_cv.out$b1.arr[,l,k],
#'                                                               X2,
#'                                                               nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_linreg_cv.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_linreg_cv.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_linreg_cv.out$eta.seq,
#'                   which.lambda.cv = grouplasso2pop_linreg_cv.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_linreg_cv.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_linreg_cv.out$which.lambda.cv.under.zero.eta,
#'                   iterations = grouplasso2pop_linreg_cv.out$iterations)
#'   
#'   class(output) <- "sempaddgt2pop_cv"
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' #' Compute semiparametric continuous-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#' #'
#' #' @param Y1 the continuous response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the continuous response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_linreg_data <- get_semipadd2pop_linreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_linreg_cv.out <- semipadd2pop_linreg_cv(Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                      X1 = semipadd2pop_linreg_data$X1,
#' #'                                                      nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                      Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                      X2 = semipadd2pop_linreg_data$X2,
#' #'                                                      nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                      w1 = 1,
#' #'                                                      w2 = 1,
#' #'                                                      w = 1,
#' #'                                                      nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                      d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                      d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                      xi = .5,
#' #'                                                      n.lambda = 5,
#' #'                                                      n.eta = 5,
#' #'                                                      lambda.min.ratio = .01,
#' #'                                                      n.folds = 5,
#' #'                                                      lambda.beta = 1,
#' #'                                                      lambda.f = 1,
#' #'                                                      eta.beta = 1,
#' #'                                                      eta.f = 1,
#' #'                                                      tol = 1e-3,
#' #'                                                      maxiter = 1000,
#' #'                                                      report.prog = FALSE)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_linreg_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_linreg_data$f1,
#' #'                                              f2 = semipadd2pop_linreg_data$f2,
#' #'                                              X1 = semipadd2pop_linreg_data$X1,
#' #'                                              X2 = semipadd2pop_linreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_linreg_cv_adapt <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,w1=1,w2=1,w=1,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_linreg_cv_adapt.out <- grouplasso2pop_linreg_cv_adapt(Y1 = Y1,
#'                                                                        X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                        groups1 = grouplasso2pop_inputs$groups1,
#'                                                                        Y2 = Y2,
#'                                                                        X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                        groups2 = grouplasso2pop_inputs$groups2,
#'                                                                        n.lambda = n.lambda,
#'                                                                        n.eta = n.eta,
#'                                                                        lambda.min.ratio = lambda.min.ratio,
#'                                                                        n.folds = n.folds,
#'                                                                        w1 = grouplasso2pop_inputs$w1,
#'                                                                        w2 = grouplasso2pop_inputs$w2,
#'                                                                        w = grouplasso2pop_inputs$w,
#'                                                                        AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                        AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                        Com = grouplasso2pop_inputs$Com,
#'                                                                        tol=tol,
#'                                                                        maxiter=maxiter,
#'                                                                        report.prog = report.prog)
#'   
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   # f1.hat.folds <- vector("list",n.lambda)
#'   # f2.hat.folds <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'     
#'     # f1.hat.folds[[l]] <- vector("list",n.eta)
#'     # f2.hat.folds[[l]] <- vector("list",n.eta)
#'     
#'     # for( k in 1:n.eta)
#'     # {
#'     #
#'     #   f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'     #   f2.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'     #
#'     # }
#'     
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1 = nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_linreg_cv_adapt.out$b1.arr[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2 = nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_linreg_cv_adapt.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#'       
#'       # for( fold in 1:n.folds)
#'       # {
#'       #
#'       #   semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'       #                                                           nonparm1 = nonparm1,
#'       #                                                           groups1 = grouplasso2pop_inputs$groups1,
#'       #                                                           knots.list1 = grouplasso2pop_inputs$knots.list1,
#'       #                                                           emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'       #                                                           QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'       #                                                           b1 = grouplasso2pop_linreg_cv_adapt.out$b1.folds.arr[,l,k,fold],
#'       #                                                           X2 = X2,
#'       #                                                           nonparm2 = nonparm2,
#'       #                                                           groups2 = grouplasso2pop_inputs$groups2,
#'       #                                                           knots.list2 = grouplasso2pop_inputs$knots.list2,
#'       #                                                           emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'       #                                                           QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'       #                                                           b2 = grouplasso2pop_linreg_cv_adapt.out$b2.folds.arr[,l,k,fold])
#'       #
#'       #   f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
#'       #   f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat
#'       #
#'       # }
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   # f1.hat.folds = f1.hat.folds,
#'                   # f2.hat.folds = f2.hat.folds,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_linreg_cv_adapt.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_linreg_cv_adapt.out$eta.seq,
#'                   lambda.initial.fit = grouplasso2pop_linreg_cv_adapt.out$lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_linreg_cv_adapt.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_linreg_cv_adapt.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_linreg_cv_adapt.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_linreg_cv_adapt.out$w1,
#'                   w2 = grouplasso2pop_linreg_cv_adapt.out$w2,
#'                   w = grouplasso2pop_linreg_cv_adapt.out$w,
#'                   iterations = grouplasso2pop_linreg_cv_adapt.out$iterations)
#'   
#'   class(output) <- "semipaddgt_cv"
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' #' Compute semiparametric continuous-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#' #'
#' #' @param Y1 the continuous response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the continuous response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.seq the sequence of lambda values
#' #' @param eta.seq the sequence of eta values
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_linreg_data <- get_semipadd2pop_linreg_data(n1 = 501, n2 = 604)
#' #'
#' #' semipadd2pop_linreg_cv_adapt.out <- semipadd2pop_linreg_cv_adapt( Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                                   X1 = semipadd2pop_linreg_data$X1,
#' #'                                                                   nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                                   Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                                   X2 = semipadd2pop_linreg_data$X2,
#' #'                                                                   nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                                   w1 = 1,
#' #'                                                                   w2 = 1,
#' #'                                                                   w = 1,
#' #'                                                                   nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                                   d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                                   d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                                   xi = .5,
#' #'                                                                   n.lambda = 5,
#' #'                                                                   n.eta = 5,
#' #'                                                                   lambda.min.ratio = .01,
#' #'                                                                   n.folds = 5,
#' #'                                                                   lambda.beta = 1,
#' #'                                                                   lambda.f = 1,
#' #'                                                                   eta.beta = 1,
#' #'                                                                   eta.f = 1,
#' #'                                                                   tol = 1e-3,
#' #'                                                                   maxiter = 1000,
#' #'                                                                   report.prog = FALSE)
#' #'
#' #' semipadd2pop_linreg_cv_adapt_fixedgrid.out <- semipadd2pop_linreg_cv_adapt_fixedgrid(Y1 = semipadd2pop_linreg_data$Y1,
#' #'                                                                                      X1 = semipadd2pop_linreg_data$X1,
#' #'                                                                                      nonparm1 = semipadd2pop_linreg_data$nonparm1,
#' #'                                                                                      Y2 = semipadd2pop_linreg_data$Y2,
#' #'                                                                                      X2 = semipadd2pop_linreg_data$X2,
#' #'                                                                                      nonparm2 = semipadd2pop_linreg_data$nonparm2,
#' #'                                                                                      nCom = semipadd2pop_linreg_data$nCom,
#' #'                                                                                      d1 = semipadd2pop_linreg_data$nonparm1*40,
#' #'                                                                                      d2 = semipadd2pop_linreg_data$nonparm2*15,
#' #'                                                                                      w1 = 1,
#' #'                                                                                      w2 = 1,
#' #'                                                                                      w = 1,
#' #'                                                                                      lambda.seq = semipadd2pop_linreg_cv_adapt.out$lambda.seq,
#' #'                                                                                      eta.seq = semipadd2pop_linreg_cv_adapt.out$eta.seq,
#' #'                                                                                      n.folds = 5,
#' #'                                                                                      lambda.beta = 1,
#' #'                                                                                      lambda.f = 1,
#' #'                                                                                      eta.beta = 1,
#' #'                                                                                      eta.f = 1,
#' #'                                                                                      tol = 1e-2,
#' #'                                                                                      maxiter = 500)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_linreg_cv_adapt_fixedgrid.out,
#' #'                        true.functions = list(f1 = semipadd2pop_linreg_data$f1,
#' #'                                              f2 = semipadd2pop_linreg_data$f2,
#' #'                                              X1 = semipadd2pop_linreg_data$X1,
#' #'                                              X2 = semipadd2pop_linreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_linreg_cv_adapt_fixedgrid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,nCom,d1,d2,xi,w1,w2,w,lambda.seq,eta.seq,lambda.initial.fit,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000)
#' {
#'   
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_linreg_cv_adapt_fixedgrid.out <- grouplasso2pop_linreg_cv_adapt_fixedgrid(Y1 = Y1,
#'                                                                                            X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                                            groups1 = grouplasso2pop_inputs$groups1,
#'                                                                                            Y2 = Y2,
#'                                                                                            X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                                            groups2 = grouplasso2pop_inputs$groups2,
#'                                                                                            lambda.seq = lambda.seq,
#'                                                                                            eta.seq = eta.seq,
#'                                                                                            lambda.initial.fit = lambda.initial.fit,
#'                                                                                            n.folds = n.folds,
#'                                                                                            w1 = grouplasso2pop_inputs$w1,
#'                                                                                            w2 = grouplasso2pop_inputs$w2,
#'                                                                                            w = grouplasso2pop_inputs$w,
#'                                                                                            AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                                            AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                                            Com = grouplasso2pop_inputs$Com,
#'                                                                                            tol = tol,
#'                                                                                            maxiter = maxiter)
#'   
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   n.lambda <- length(lambda.seq)
#'   n.eta <- length(eta.seq)
#'   
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'     
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$b1.arr[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = lambda.seq,
#'                   eta.seq = eta.seq,
#'                   lambda.initial.fit = lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$w1,
#'                   w2 = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$w2,
#'                   w = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$w,
#'                   iterations = grouplasso2pop_linreg_cv_adapt_fixedgrid.out$iterations)
#'   
#'   class(output) <- "semipaddgt_cv"
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 the binary response vector of data set 1
#' #' @param XX1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the binary response vector of data set 2
#' #' @param XX2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_logreg.out <- semipadd2pop_logreg(Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                X1 = semipadd2pop_logreg_data$X1,
#' #'                                                nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                X2 = semipadd2pop_logreg_data$X2,
#' #'                                                nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                rho1 = 2,
#' #'                                                rho2 = 1,
#' #'                                                w1 = 1,
#' #'                                                w2 = 1,
#' #'                                                w = 1,
#' #'                                                nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                d1 = semipadd2pop_logreg_data$nonparm1*25,
#' #'                                                d2 = semipadd2pop_logreg_data$nonparm2*15,
#' #'                                                xi = .5,
#' #'                                                lambda.beta = .01,
#' #'                                                lambda.f = .01,
#' #'                                                eta.beta = .01,
#' #'                                                eta.f = .01,
#' #'                                                tol = 1e-3,
#' #'                                                maxiter = 500,
#' #'                                                plot_obj = FALSE)
#' #'
#' #' plot_semipaddgt2pop(semipadd2pop_logreg.out,
#' #'                     true.functions=list(f1 = semipadd2pop_logreg_data$f1,
#' #'                                         f2 = semipadd2pop_logreg_data$f2,
#' #'                                         X1 = semipadd2pop_logreg_data$X1,
#' #'                                         X2 = semipadd2pop_logreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_logreg <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,tol=1e-4,maxiter=500,plot_obj=FALSE)
#' {
#' 
#'   # prepare input for grouplasso2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#' 
#'   # get group lasso estimators
#'   grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
#'                                                      rX1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                      groups1 = grouplasso2pop_inputs$groups1,
#'                                                      rY2 = Y2,
#'                                                      rX2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                      groups2 = grouplasso2pop_inputs$groups2,
#'                                                      rho1 = rho1,
#'                                                      rho2 = rho2,
#'                                                      lambda = grouplasso2pop_inputs$lambda,
#'                                                      eta = grouplasso2pop_inputs$eta,
#'                                                      w1 = grouplasso2pop_inputs$w1,
#'                                                      w2 = grouplasso2pop_inputs$w2,
#'                                                      w = grouplasso2pop_inputs$w,
#'                                                      rAA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                      rAA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                      rCom = grouplasso2pop_inputs$Com,
#'                                                      tol = tol,
#'                                                      maxiter = maxiter)
#' 
#'   # construct fitted functions from grouplasso2pop output
#'   semipadd2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                         nonparm1 = nonparm1,
#'                                                         groups1 = grouplasso2pop_inputs$groups1,
#'                                                         knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                         emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                         QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                         b1 = grouplasso2pop_logreg.out$beta1.hat,
#'                                                         X2 = X2,
#'                                                         nonparm2 = nonparm2,
#'                                                         groups2 = grouplasso2pop_inputs$groups2,
#'                                                         knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                         emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                         QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                         b2 = grouplasso2pop_logreg.out$beta2.hat)
#' 
#'   # collect output
#'   output <- list( f1.hat = semipadd2pop_fitted$f1.hat,
#'                   f2.hat = semipadd2pop_fitted$f2.hat,
#'                   f1.hat.design = semipadd2pop_fitted$f1.hat.design,
#'                   f2.hat.design = semipadd2pop_fitted$f2.hat.design,
#'                   beta1.hat = semipadd2pop_fitted$beta1.hat,
#'                   beta2.hat = semipadd2pop_fitted$beta2.hat,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com)
#' 
#'   class(output) <- "semipaddgt2pop"
#' 
#'   return(output)
#' }
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity over a grid of tuning parameter values
#' #'
#' #' @param Y1 the binary response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the binary response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_logreg_grid.out <- semipadd2pop_logreg_grid(Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                          X1 = semipadd2pop_logreg_data$X1,
#' #'                                                          nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                          Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                          X2 = semipadd2pop_logreg_data$X2,
#' #'                                                          nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                          rho1 = 2,
#' #'                                                          rho2 = 1, 
#' #'                                                          nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                          d1 = semipadd2pop_logreg_data$nonparm1*25,
#' #'                                                          d2 = semipadd2pop_logreg_data$nonparm2*15,
#' #'                                                          xi = .5,
#' #'                                                          w1 = 1,
#' #'                                                          w2 = 1,
#' #'                                                          w = 1,
#' #'                                                          lambda.beta = 1,
#' #'                                                          lambda.f = 1,
#' #'                                                          eta.beta = 1,
#' #'                                                          eta.f = 1,
#' #'                                                          n.lambda = 5,
#' #'                                                          n.eta = 5,
#' #'                                                          lambda.min.ratio = .01,
#' #'                                                          tol = 1e-3,
#' #'                                                          maxiter = 500,
#' #'                                                          report.prog = TRUE)
#' #'
#' #' plot_semipaddgt2pop_grid(semipadd2pop_logreg_grid.out,
#' #'                          true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#' #'                                                f2 = semipadd2pop_logreg_data$f2,
#' #'                                                X1 = semipadd2pop_logreg_data$X1,
#' #'                                                X2 = semipadd2pop_logreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_logreg_grid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#' 
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#' 
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = Y1,
#'                                                                X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                Y2 = Y2,
#'                                                                X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                rho1 = rho1,
#'                                                                rho2 = rho2,
#'                                                                n.lambda = n.lambda,
#'                                                                n.eta = n.eta,
#'                                                                lambda.min.ratio = lambda.min.ratio,
#'                                                                lambda.max.ratio = lambda.max.ratio,
#'                                                                w1 = grouplasso2pop_inputs$w1,
#'                                                                w2 = grouplasso2pop_inputs$w2,
#'                                                                w = grouplasso2pop_inputs$w,
#'                                                                AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                Com = grouplasso2pop_inputs$Com,
#'                                                                tol = tol,
#'                                                                maxiter = maxiter,
#'                                                                report.prog = report.prog)
#' 
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#' 
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#' 
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'   }
#' 
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#' 
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1 = nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_logreg_grid.out$b1[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2 = nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_logreg_grid.out$b2[,l,k])
#' 
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#' 
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#' 
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#' 
#'     }
#' 
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   P1.hat = grouplasso2pop_logreg_grid.out$P1.arr,
#'                   P2.hat = grouplasso2pop_logreg_grid.out$P2.arr,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   lambda.seq = grouplasso2pop_logreg_grid.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_logreg_grid.out$eta.seq,
#'                   iterations = grouplasso2pop_logreg_grid.out$iterations)
#' 
#'   class(output) <- "semipaddgt2pop_grid"
#' 
#'   return(output)
#' 
#' }
#' 
#' 
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters
#' #'
#' #' @param Y1 the binary response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the binary response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_logreg_cv.out <- semipadd2pop_logreg_cv(Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                      X1 = semipadd2pop_logreg_data$X1,
#' #'                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                      Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                      X2 = semipadd2pop_logreg_data$X2,
#' #'                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                      rho1 = 2,
#' #'                                                      rho2 = 1,
#' #'                                                      w1 = 1,
#' #'                                                      w2 = 1,
#' #'                                                      w = 1,
#' #'                                                      nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                      d1 = semipadd2pop_logreg_data$nonparm1*25,
#' #'                                                      d2 = semipadd2pop_logreg_data$nonparm2*15,
#' #'                                                      xi = .5,
#' #'                                                      n.lambda = 5,
#' #'                                                      n.eta = 5,
#' #'                                                      lambda.min.ratio = .01,
#' #'                                                      n.folds = 5,
#' #'                                                      lambda.beta = 1,
#' #'                                                      lambda.f = 1,
#' #'                                                      eta.beta = 1,
#' #'                                                      eta.f = 1,
#' #'                                                      tol = 1e-3,
#' #'                                                      maxiter = 1000,
#' #'                                                      report.prog = FALSE)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#' #'                                              f2 = semipadd2pop_logreg_data$f2,
#' #'                                              X1 = semipadd2pop_logreg_data$X1,
#' #'                                              X2 = semipadd2pop_logreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_logreg_cv <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#' 
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#' 
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_logreg_cv.out <- grouplasso2pop_logreg_cv(Y1 = Y1,
#'                                                            X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                            groups1 = grouplasso2pop_inputs$groups1,
#'                                                            Y2 = Y2,
#'                                                            X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                            groups2 = grouplasso2pop_inputs$groups2,
#'                                                            rho1 = rho1,
#'                                                            rho2 = rho2,
#'                                                            n.lambda = n.lambda,
#'                                                            n.eta = n.eta,
#'                                                            lambda.min.ratio = lambda.min.ratio,
#'                                                            lambda.max.ratio = lambda.max.ratio,
#'                                                            n.folds = n.folds,
#'                                                            w1 = grouplasso2pop_inputs$w1,
#'                                                            w2 = grouplasso2pop_inputs$w2,
#'                                                            w = grouplasso2pop_inputs$w,
#'                                                            AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                            AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                            Com = grouplasso2pop_inputs$Com,
#'                                                            tol = tol,
#'                                                            maxiter = maxiter,
#'                                                            report.prog = report.prog)
#' 
#' 
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#' 
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#' 
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#' 
#'   f1.hat.folds <- vector("list",n.lambda)
#'   f2.hat.folds <- vector("list",n.lambda)
#' 
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#' 
#'     f1.hat.folds[[l]] <- vector("list",n.eta)
#'     f2.hat.folds[[l]] <- vector("list",n.eta)
#' 
#'     for( k in 1:n.eta)
#'     {
#' 
#'       f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       f2.hat.folds[[l]][[k]] <- vector("list",n.folds)
#' 
#'     }
#' 
#'   }
#' 
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#' 
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1 = nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_logreg_cv.out$b1.arr[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2 = nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_logreg_cv.out$b2.arr[,l,k])
#' 
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#' 
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#' 
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#' 
#'       for( fold in 1:n.folds)
#'       {
#' 
#'         semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                 nonparm1 = nonparm1,
#'                                                                 groups1 = grouplasso2pop_inputs$groups1,
#'                                                                 knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                 emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                 QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                 b1 = grouplasso2pop_logreg_cv.out$b1.folds.arr[,l,k,fold],
#'                                                                 X2 = X2,
#'                                                                 nonparm2 = nonparm2,
#'                                                                 groups2 = grouplasso2pop_inputs$groups2,
#'                                                                 knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                 emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                 QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                 b2 = grouplasso2pop_logreg_cv.out$b2.folds.arr[,l,k,fold])
#' 
#'         f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
#'         f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat
#' 
#'       }
#' 
#'     }
#' 
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.folds = f1.hat.folds,
#'                   f2.hat.folds = f2.hat.folds,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   rho1 = rho1, 
#'                   rho2 = rho2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_logreg_cv.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_logreg_cv.out$eta.seq,
#'                   which.lambda.cv = grouplasso2pop_logreg_cv.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_logreg_cv.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv.out$which.lambda.cv.under.zero.eta,
#'                   iterations = grouplasso2pop_logreg_cv.out$iterations)
#' 
#'   class(output) <- "sempaddgt2pop_cv"
#' 
#'   return(output)
#' 
#' }
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#' #'
#' #' @param Y1 the binary response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the binary response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#' #'                                                          n2 = 604)
#' #'
#' #' semipadd2pop_logreg_cv.out <- semipadd2pop_logreg_cv(Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                      X1 = semipadd2pop_logreg_data$X1,
#' #'                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                      Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                      X2 = semipadd2pop_logreg_data$X2,
#' #'                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                      rho1 = 2,
#' #'                                                      rho2 = 1,
#' #'                                                      w1 = 1,
#' #'                                                      w2 = 1,
#' #'                                                      w = 1,
#' #'                                                      nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                      d1 = semipadd2pop_logreg_data$nonparm1*25,
#' #'                                                      d2 = semipadd2pop_logreg_data$nonparm2*15,
#' #'                                                      xi = .5,
#' #'                                                      n.lambda = 5,
#' #'                                                      n.eta = 5,
#' #'                                                      lambda.min.ratio = .001,
#' #'                                                      n.folds = 5,
#' #'                                                      lambda.beta = 1,
#' #'                                                      lambda.f = 1,
#' #'                                                      eta.beta = 1,
#' #'                                                      eta.f = 1,
#' #'                                                      tol = 1e-3,
#' #'                                                      maxiter = 1000,
#' #'                                                      report.prog = FALSE)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#' #'                                              f2 = semipadd2pop_logreg_data$f2,
#' #'                                              X1 = semipadd2pop_logreg_data$X1,
#' #'                                              X2 = semipadd2pop_logreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_logreg_cv_adapt <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1=1,w2=1,w=1,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio = .01,lambda.max.ratio=1,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#' 
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#' 
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_logreg_cv_adapt.out <- grouplasso2pop_logreg_cv_adapt(Y1 = Y1,
#'                                                                        X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                        groups1 = grouplasso2pop_inputs$groups1,
#'                                                                        Y2 = Y2,
#'                                                                        X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                        groups2 = grouplasso2pop_inputs$groups2,
#'                                                                        rho1 = rho1,
#'                                                                        rho2 = rho2,
#'                                                                        n.lambda = n.lambda,
#'                                                                        n.eta = n.eta,
#'                                                                        lambda.min.ratio = lambda.min.ratio,
#'                                                                        lambda.max.ratio = lambda.max.ratio,
#'                                                                        n.folds = n.folds,
#'                                                                        w1 = grouplasso2pop_inputs$w1,
#'                                                                        w2 = grouplasso2pop_inputs$w2,
#'                                                                        w = grouplasso2pop_inputs$w,
#'                                                                        AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                        AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                        Com = grouplasso2pop_inputs$Com,
#'                                                                        tol = tol,
#'                                                                        maxiter = maxiter,
#'                                                                        report.prog = report.prog)
#' 
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#' 
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#' 
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#' 
#'   f1.hat.folds <- vector("list",n.lambda)
#'   f2.hat.folds <- vector("list",n.lambda)
#' 
#'   for(l in 1:n.lambda)
#'   {
#'     
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#' 
#'     f1.hat.folds[[l]] <- vector("list",n.eta)
#'     f2.hat.folds[[l]] <- vector("list",n.eta)
#' 
#'     for( k in 1:n.eta)
#'     {
#' 
#'       f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       f2.hat.folds[[l]][[k]] <- vector("list",n.folds)
#' 
#'     }
#' 
#'   }
#' 
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#' 
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1 = nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_logreg_cv_adapt.out$b1.arr[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2 = nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_logreg_cv_adapt.out$b2.arr[,l,k])
#' 
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#' 
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#' 
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#' 
#'       # for( fold in 1:n.folds)
#'       # {
#'       # 
#'       #   semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'       #                                                           nonparm1 = nonparm1,
#'       #                                                           groups1 = grouplasso2pop_inputs$groups1,
#'       #                                                           knots.list1 = grouplasso2pop_inputs$knots.list1,
#'       #                                                           emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'       #                                                           QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'       #                                                           b1 = grouplasso2pop_logreg_cv_adapt.out$b1.folds.arr[,l,k,fold],
#'       #                                                           X2 = X2,
#'       #                                                           nonparm2 = nonparm2,
#'       #                                                           groups2 = grouplasso2pop_inputs$groups2,
#'       #                                                           knots.list2 = grouplasso2pop_inputs$knots.list2,
#'       #                                                           emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'       #                                                           QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'       #                                                           b2 = grouplasso2pop_logreg_cv_adapt.out$b2.folds.arr[,l,k,fold])
#'       # 
#'       #   f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
#'       #   f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat
#'       # 
#'       # }
#' 
#'     }
#' 
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.folds = f1.hat.folds,
#'                   f2.hat.folds = f2.hat.folds,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   P1.hat = grouplasso2pop_logreg_cv_adapt.out$P1.arr,
#'                   P2.hat = grouplasso2pop_logreg_cv_adapt.out$P2.arr,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_logreg_cv_adapt.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_logreg_cv_adapt.out$eta.seq,
#'                   lambda.initial.fit = grouplasso2pop_logreg_cv_adapt.out$lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_logreg_cv_adapt.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_logreg_cv_adapt.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_adapt.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_logreg_cv_adapt.out$w1,
#'                   w2 = grouplasso2pop_logreg_cv_adapt.out$w2,
#'                   w = grouplasso2pop_logreg_cv_adapt.out$w,
#'                   iterations = grouplasso2pop_logreg_cv_adapt.out$iterations)
#' 
#'   class(output) <- "semipaddgt_cv"
#' 
#'   return(output)
#' 
#' }
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#' #'
#' #' @param Y1 the binary response vector of data set 1
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 the binary response vector of data set 2
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{XX1} and \code{XX2} after the column of ones corresponding to the intercept.
#' #' @param d the dimension of the B-spline basis to be used when fitting the nonparametric effects
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.seq the sequence of lambda values
#' #' @param eta.seq the sequence of eta values
#' #' @param n.folds the number of crossvalidation folds
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations
#' #' @return Returns the estimator of the semiparametric additive model
#' #'
#' #' @examples
#' #' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501, n2 = 604)
#' #'
#' #' semipadd2pop_logreg_cv_adapt.out <- semipadd2pop_logreg_cv_adapt( Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                                   X1 = semipadd2pop_logreg_data$X1,
#' #'                                                                   nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                                   Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                                   X2 = semipadd2pop_logreg_data$X2,
#' #'                                                                   nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                                   w1 = 1,
#' #'                                                                   w2 = 1,
#' #'                                                                   w = 1,
#' #'                                                                   nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                                   d1 = semipadd2pop_logreg_data$nonparm1*25,
#' #'                                                                   d2 = semipadd2pop_logreg_data$nonparm2*15,
#' #'                                                                   xi = .5,
#' #'                                                                   n.lambda = 5,
#' #'                                                                   n.eta = 5,
#' #'                                                                   lambda.min.ratio = .01,
#' #'                                                                   n.folds = 5,
#' #'                                                                   lambda.beta = 1,
#' #'                                                                   lambda.f = 1,
#' #'                                                                   eta.beta = 1,
#' #'                                                                   eta.f = 1,
#' #'                                                                   tol = 1e-3,
#' #'                                                                   maxiter = 1000,
#' #'                                                                   report.prog = FALSE)
#' #'
#' #' semipadd2pop_logreg_cv_adapt_fixedgrid.out <- semipadd2pop_logreg_cv_adapt_fixedgrid(Y1 = semipadd2pop_logreg_data$Y1,
#' #'                                                                                      X1 = semipadd2pop_logreg_data$X1,
#' #'                                                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#' #'                                                                                      Y2 = semipadd2pop_logreg_data$Y2,
#' #'                                                                                      X2 = semipadd2pop_logreg_data$X2,
#' #'                                                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#' #'                                                                                      rho1 = 2,
#' #'                                                                                      rho2 = 1,
#' #'                                                                                      nCom = semipadd2pop_logreg_data$nCom,
#' #'                                                                                      d = 20,
#' #'                                                                                      xi = .5,
#' #'                                                                                      w1 = 1,
#' #'                                                                                      w2 = 1,
#' #'                                                                                      w = 1,
#' #'                                                                                      lambda.seq = semipadd2pop_logreg_cv_adapt.out$lambda.seq,
#' #'                                                                                      eta.seq = semipadd2pop_logreg_cv_adapt.out$eta.seq,
#' #'                                                                                      n.folds = 5,
#' #'                                                                                      lambda.beta = 1,
#' #'                                                                                      lambda.f = 1,
#' #'                                                                                      eta.beta = 1,
#' #'                                                                                      eta.f = 1,
#' #'                                                                                      tol = 1e-2,
#' #'                                                                                      maxiter = 500)
#' #'
#' #' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv_adapt_fixedgrid.out,
#' #'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#' #'                                              f2 = semipadd2pop_logreg_data$f2,
#' #'                                              X1 = semipadd2pop_logreg_data$X1,
#' #'                                              X2 = semipadd2pop_logreg_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_logreg_cv_adapt_fixedgrid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,nCom,d1,d2,xi,w1,w2,w,lambda.seq,eta.seq,lambda.initial.fit,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#' 
#'   # prepare input for grouplassogt2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#' 
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_logreg_cv_adapt_fixedgrid.out <- grouplasso2pop_logreg_cv_adapt_fixedgrid(Y1 = Y1,
#'                                                                                            X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                                            groups1 = grouplasso2pop_inputs$groups1,
#'                                                                                            Y2 = Y2,
#'                                                                                            X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                                            groups2 = grouplasso2pop_inputs$groups2,
#'                                                                                            rho1 = rho1,
#'                                                                                            rho2 = rho2,
#'                                                                                            lambda.seq = lambda.seq,
#'                                                                                            eta.seq = eta.seq,
#'                                                                                            lambda.initial.fit = lambda.initial.fit,
#'                                                                                            n.folds = n.folds,
#'                                                                                            w1 = grouplasso2pop_inputs$w1,
#'                                                                                            w2 = grouplasso2pop_inputs$w2,
#'                                                                                            w = grouplasso2pop_inputs$w,
#'                                                                                            AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                                            AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                                            Com = grouplasso2pop_inputs$Com,
#'                                                                                            tol = tol,
#'                                                                                            maxiter = maxiter,
#'                                                                                            report.prog = report.prog)
#'   
#'   which.lambda.cv <- grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv
#'   which.eta.cv <- grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.eta.cv
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#' 
#'   n.lambda <- length(lambda.seq)
#'   n.eta <- length(eta.seq)
#' 
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#' 
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#' 
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#' 
#'   f1.hat.folds <- vector("list",n.lambda)
#'   f2.hat.folds <- vector("list",n.lambda)
#' 
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#' 
#'   }
#' 
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#' 
#'       semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                               nonparm1 = nonparm1,
#'                                                               groups1 = grouplasso2pop_inputs$groups1,
#'                                                               knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                               emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                               QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                               b1 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$b1.arr[,l,k],
#'                                                               X2 = X2,
#'                                                               nonparm2,
#'                                                               groups2 = grouplasso2pop_inputs$groups2,
#'                                                               knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                               emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                               QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                               b2 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$b2.arr[,l,k])
#' 
#'       f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat
#' 
#'       f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design
#' 
#'       beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat
#' 
#'     }
#' 
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   P1.hat = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$P1.arr,
#'                   P2.hat = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$P2.arr,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = lambda.seq,
#'                   eta.seq = eta.seq,
#'                   lambda.initial.fit = lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w1,
#'                   w2 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w2,
#'                   w = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w,
#'                   iterations = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$iterations)
#' 
#'   class(output) <- "semipadd2pop_gt_cv"
#' 
#'   return(output)
#' 
#' }
#' 
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se1 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp1 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se2 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp2 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common
#' #' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations (EM steps)
#' #' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' #' @return Returns the estimator of the semiparametric additive model with group testing data
#' #'
#' #' @examples
#' #' semipadd2pop_gt_data <- get_semipadd2pop_gt_data(n1 = 500, n2 = 604)
#' #'
#' #' semipadd2pop_gt.out <- semipadd2pop_gt(Y1 = semipadd2pop_gt_data$Y1,
#' #'                                      Z1 = semipadd2pop_gt_data$Z1,
#' #'                                      Se1 = semipadd2pop_gt_data$Se1,
#' #'                                      Sp1 = semipadd2pop_gt_data$Sp1,
#' #'                                      X1 = semipadd2pop_gt_data$X1,
#' #'                                      nonparm1 = semipadd2pop_gt_data$nonparm1,
#' #'                                      Y2 = semipadd2pop_gt_data$Y2,
#' #'                                      Z2 = semipadd2pop_gt_data$Z2,
#' #'                                      Se2 = semipadd2pop_gt_data$Se2,
#' #'                                      Sp2 = semipadd2pop_gt_data$Sp2,
#' #'                                      X2 = semipadd2pop_gt_data$X2,
#' #'                                      nonparm2 = semipadd2pop_gt_data$nonparm2,
#' #'                                      rho1 = 2,
#' #'                                      rho2 = 1,
#' #'                                      w1 = 1,
#' #'                                      w2 = 1,
#' #'                                      w = 1,
#' #'                                      nCom = 4,
#' #'                                      d1 = semipadd2pop_gt_data$nonparm1 * 15,
#' #'                                      d2 = semipadd2pop_gt_data$nonparm2 * 10,
#' #'                                      xi = 1,
#' #'                                      lambda.beta = 1,
#' #'                                      lambda.f = 1,
#' #'                                      eta.beta = 1,
#' #'                                      eta.f = 1,
#' #'                                      tol = 1e-2,
#' #'                                      maxiter = 500,
#' #'                                      report.prog = TRUE)
#' #'
#' #' plot_semipadd2pop_gt(semipadd2pop_gt.out,
#' #'                     true.functions=list(f1 = semipadd2pop_gt_data$f1,
#' #'                                         f2 = semipadd2pop_gt_data$f2,
#' #'                                         X1 = semipadd2pop_gt_data$X1,
#' #'                                         X2 = semipadd2pop_gt_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_gt <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,E.approx = FALSE,tol=1e-3,maxiter=500,report.prog=FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop_gt function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimates using the EM-algorithm
#'   grouplasso2pop_gt.out <- grouplasso2pop_gt(Y1 = Y1,
#'                                              Z1 = Z1,
#'                                              Se1 = Se1,
#'                                              Sp1 = Sp1,
#'                                              X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                              groups1 = grouplasso2pop_inputs$groups1,
#'                                              Y2 = Y2,
#'                                              Z2 = Z2,
#'                                              Se2 = Se2,
#'                                              Sp2 = Sp2,
#'                                              X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                              groups2 = grouplasso2pop_inputs$groups2,
#'                                              rho1 = rho1,
#'                                              rho2 = rho2,
#'                                              lambda = grouplasso2pop_inputs$lambda,
#'                                              eta = grouplasso2pop_inputs$eta,
#'                                              w1 = grouplasso2pop_inputs$w1,
#'                                              w2 = grouplasso2pop_inputs$w2,
#'                                              w = grouplasso2pop_inputs$w,
#'                                              AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                              AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                              Com = grouplasso2pop_inputs$Com,
#'                                              E.approx = E.approx,
#'                                              tol = tol,
#'                                              maxiter = maxiter,
#'                                              init = NULL,
#'                                              report.prog = report.prog)
#'   
#'   # construct fitted functions from grouplasso2pop output
#'   semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                            nonparm1 = nonparm1,
#'                                                            groups1 = grouplasso2pop_inputs$groups1,
#'                                                            knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                            emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                            QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                            b1 = grouplasso2pop_gt.out$beta1.hat,
#'                                                            X2 = X2,
#'                                                            nonparm2 = nonparm2,
#'                                                            groups2 = grouplasso2pop_inputs$groups2,
#'                                                            knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                            emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                            QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                            b2 = grouplasso2pop_gt.out$beta2.hat)
#'   
#'   
#'   # collect output
#'   output <- list( f1.hat = semipadd2pop_gt_fitted$f1.hat,
#'                   f2.hat = semipadd2pop_gt_fitted$f2.hat,
#'                   f1.hat.design = semipadd2pop_gt_fitted$f1.hat.design,
#'                   f2.hat.design = semipadd2pop_gt_fitted$f2.hat.design,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   inner.iter = grouplasso2pop_gt.out$inner.iter)
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' 
#' #' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se1 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp1 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se2 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp2 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the  B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations (EM steps)
#' #' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' #' @return Returns the estimator of the semiparametric additive model with group testing data
#' #'
#' #' @examples
#' #' semipadd2pop_gt_data <- get_semipadd2pop_gt_data(n1 = 500, n2 = 604)
#' #'
#' #' semipadd2pop_gt_grid.out <- semipadd2pop_gt_grid(Y1 = semipadd2pop_gt_data$Y1,
#' #'                                                Z1 = semipadd2pop_gt_data$Z1,
#' #'                                                Se1 = semipadd2pop_gt_data$Se1,
#' #'                                                Sp1 = semipadd2pop_gt_data$Sp1,
#' #'                                                X1 = semipadd2pop_gt_data$X1,
#' #'                                                nonparm1 = semipadd2pop_gt_data$nonparm1,
#' #'                                                Y2 = semipadd2pop_gt_data$Y2,
#' #'                                                Z2 = semipadd2pop_gt_data$Z2,
#' #'                                                Se2 = semipadd2pop_gt_data$Se2,
#' #'                                                Sp2 = semipadd2pop_gt_data$Sp2,
#' #'                                                X2 = semipadd2pop_gt_data$X2,
#' #'                                                nonparm2 = semipadd2pop_gt_data$nonparm2,
#' #'                                                rho1 = 2,
#' #'                                                rho2 = 1,
#' #'                                                w1 = 1,
#' #'                                                w2 = 1,
#' #'                                                w = 1,
#' #'                                                nCom = 4,
#' #'                                                d1 = semipadd2pop_gt_data$nonparm1 * 15,
#' #'                                                d2 = semipadd2pop_gt_data$nonparm2 * 10,
#' #'                                                xi = 1,
#' #'                                                n.lambda = 3,
#' #'                                                n.eta = 3,
#' #'                                                lambda.min.ratio=.01,
#' #'                                                lambda.beta = 1,
#' #'                                                lambda.f = 1,
#' #'                                                eta.beta = 1,
#' #'                                                eta.f = 1,
#' #'                                                tol = 1e-2,
#' #'                                                maxiter = 500,
#' #'                                                report.prog = TRUE)
#' #'
#' #' plot_semipadd2pop_gt_grid(semipadd2pop_gt_grid.out,
#' #'                          true.functions = list(f1 = semipadd2pop_gt_data$f1,
#' #'                                                f2 = semipadd2pop_gt_data$f2,
#' #'                                                X1 = semipadd2pop_gt_data$X1,
#' #'                                                X2 = semipadd2pop_gt_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_gt_grid <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w2,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values
#'   grouplasso2pop_gt_grid.out <- grouplasso2pop_gt_grid(Y1 = Y1,
#'                                                        Z1 = Z1,
#'                                                        Se1 = Se1,
#'                                                        Sp1 = Sp1,
#'                                                        X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                        groups1 = grouplasso2pop_inputs$groups1,
#'                                                        Y2 = Y2,
#'                                                        Z2 = Z2,
#'                                                        Se2 = Se2,
#'                                                        Sp2 = Sp2,
#'                                                        X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                        groups2 =  grouplasso2pop_inputs$groups2,
#'                                                        rho1 = rho1,
#'                                                        rho2 = rho2,
#'                                                        n.lambda = n.lambda,
#'                                                        n.eta = n.eta,
#'                                                        lambda.min.ratio = lambda.min.ratio,
#'                                                        lambda.max.ratio = lambda.max.ratio,
#'                                                        w1 = grouplasso2pop_inputs$w1,
#'                                                        w2 = grouplasso2pop_inputs$w2,
#'                                                        w = grouplasso2pop_inputs$w,
#'                                                        AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                        AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                        Com = grouplasso2pop_inputs$Com,
#'                                                        E.approx = E.approx,
#'                                                        tol = tol,
#'                                                        maxiter = maxiter,
#'                                                        report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                nonparm1 = nonparm1,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                b1 = grouplasso2pop_gt_grid.out$b1.arr[,l,k],
#'                                                                X2 = X2,
#'                                                                nonparm2 = nonparm2,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                b2 = grouplasso2pop_gt_grid.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f2.hat.design
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   P1.hat = grouplasso2pop_gt_grid.out$P1.arr,
#'                   P2.hat = grouplasso2pop_gt_grid.out$P2.arr,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   lambda.seq = grouplasso2pop_gt_grid.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_gt_grid.out$eta.seq,
#'                   iterations = grouplasso2pop_gt_grid.out$iterations)
#'   
#'   return(output)
#'   
#' }
#' 
#' #' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se1 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp1 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se2 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp2 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations (EM steps)
#' #' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' #' @return Returns the estimator of the semiparametric additive model with group testing data
#' #'
#' #' @examples
#' #' semipadd2pop_gt_data <- get_semipadd2pop_gt_data(n1 = 500, n2 = 604)
#' #'
#' #' semipadd2pop_gt_cv.out <- semipadd2pop_gt_cv(Y1 = semipadd2pop_gt_data$Y1,
#' #'                                            Z1 = semipadd2pop_gt_data$Z1,
#' #'                                            Se1 = semipadd2pop_gt_data$Se1,
#' #'                                            Sp1 = semipadd2pop_gt_data$Sp1,
#' #'                                            X1 = semipadd2pop_gt_data$X1,
#' #'                                            nonparm1 = semipadd2pop_gt_data$nonparm1,
#' #'                                            Y2 = semipadd2pop_gt_data$Y2,
#' #'                                            Z2 = semipadd2pop_gt_data$Z2,
#' #'                                            Se2 = semipadd2pop_gt_data$Se2,
#' #'                                            Sp2 = semipadd2pop_gt_data$Sp2,
#' #'                                            X2 = semipadd2pop_gt_data$X2,
#' #'                                            nonparm2 = semipadd2pop_gt_data$nonparm2,
#' #'                                            rho1 = 2,
#' #'                                            rho2 = 1,
#' #'                                            w1 = 1,
#' #'                                            w2 = 1,
#' #'                                            w = 1,
#' #'                                            nCom = 4,
#' #'                                            d1 = semipadd2pop_gt_data$nonparm1 * 15,
#' #'                                            d2 = semipadd2pop_gt_data$nonparm2 * 10,
#' #'                                            xi = 1,
#' #'                                            n.lambda = 3,
#' #'                                            n.eta = 3,
#' #'                                            lambda.min.ratio =.001,
#' #'                                            lambda.beta = 1,
#' #'                                            lambda.f = 1,
#' #'                                            eta.beta = 1,
#' #'                                            eta.f = 1,
#' #'                                            tol = 1e-2,
#' #'                                            maxiter = 500,
#' #'                                            report.prog = TRUE)
#' #'
#' #' plot_semipadd2pop_gt_cv(semipadd2pop_gt_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_gt_data$f1,
#' #'                                              f2 = semipadd2pop_gt_data$f2,
#' #'                                              X1 = semipadd2pop_gt_data$X1,
#' #'                                              X2 = semipadd2pop_gt_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_gt_cv <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop_gt function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
#'   grouplasso2pop_gt_cv.out <- grouplasso2pop_gt_cv(Y1 = Y1,
#'                                                    Z1 = Z1,
#'                                                    Se1 = Se1,
#'                                                    Sp1 = Sp1,
#'                                                    X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                    groups1 = grouplasso2pop_inputs$groups1,
#'                                                    Y2 = Y2,
#'                                                    Z2 = Z2,
#'                                                    Se2 = Se2,
#'                                                    Sp2 = Sp2,
#'                                                    X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                    groups2 = grouplasso2pop_inputs$groups2,
#'                                                    rho1 = rho1,
#'                                                    rho2 = rho2,
#'                                                    n.lambda = n.lambda,
#'                                                    n.eta = n.eta,
#'                                                    lambda.min.ratio = lambda.min.ratio,
#'                                                    lambda.max.ratio = lambda.max.ratio,
#'                                                    n.folds = n.folds,
#'                                                    w1 = grouplasso2pop_inputs$w1,
#'                                                    w2 = grouplasso2pop_inputs$w2,
#'                                                    w = grouplasso2pop_inputs$w,
#'                                                    AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                    AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                    Com = grouplasso2pop_inputs$Com,
#'                                                    E.approx = E.approx,
#'                                                    tol = tol,
#'                                                    maxiter = maxiter,
#'                                                    report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   f1.hat.folds <- vector("list",n.lambda)
#'   f2.hat.folds <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'     
#'     f1.hat.folds[[l]] <- vector("list",n.eta)
#'     f2.hat.folds[[l]] <- vector("list",n.eta)
#'     
#'     for( k in 1:n.eta)
#'     {
#'       
#'       f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       f2.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       
#'     }
#'     
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                nonparm1 = nonparm1,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                b1 = grouplasso2pop_gt_cv.out$b1.arr[,l,k],
#'                                                                X2 = X2,
#'                                                                nonparm2 = nonparm2,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                b2 = grouplasso2pop_gt_cv.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipadd2pop_gt_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipadd2pop_gt_fitted$beta2.hat
#'       
#'       for( fold in 1:n.folds)
#'       {
#'         
#'         semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                  nonparm1 = nonparm1,
#'                                                                  groups1 = grouplasso2pop_inputs$groups1,
#'                                                                  knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                  emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                  QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                  b1 = grouplasso2pop_gt_cv.out$b1.folds.arr[,l,k,fold],
#'                                                                  X2 = X2,
#'                                                                  nonparm2 = nonparm2,
#'                                                                  groups2 = grouplasso2pop_inputs$groups2,
#'                                                                  knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                  emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                  QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                  b2 = grouplasso2pop_gt_cv.out$b2.folds.arr[,l,k,fold])
#'         
#'         f1.hat.folds[[l]][[k]][[fold]] <- semipadd2pop_gt_fitted$f1.hat
#'         f2.hat.folds[[l]][[k]][[fold]] <- semipadd2pop_gt_fitted$f2.hat
#'         
#'       }
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.folds = f1.hat.folds,
#'                   f2.hat.folds = f2.hat.folds,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_gt_cv.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_gt_cv.out$eta.seq,
#'                   which.lambda.cv = grouplasso2pop_gt_cv.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_gt_cv.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_gt_cv.out$which.lambda.cv.under.zero.eta,
#'                   iterations = grouplasso2pop_gt_cv.out$iterations)
#'   
#'   class(output) <- "semipadd2pop_gt_cv"
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' 
#' #' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se1 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp1 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se2 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp2 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations (EM steps)
#' #' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' #' @return Returns the estimator of the semiparametric additive model with group testing data
#' #'
#' #' @examples
#' #' semipadd2pop_gt_data <- get_semipadd2pop_gt_data(n1 = 500, n2 = 604)
#' #'
#' #' semipadd2pop_gt_cv.out <- semipadd2pop_gt_cv(Y1 = semipadd2pop_gt_data$Y1,
#' #'                                            Z1 = semipadd2pop_gt_data$Z1,
#' #'                                            Se1 = semipadd2pop_gt_data$Se1,
#' #'                                            Sp1 = semipadd2pop_gt_data$Sp1,
#' #'                                            X1 = semipadd2pop_gt_data$X1,
#' #'                                            nonparm1 = semipadd2pop_gt_data$nonparm1,
#' #'                                            Y2 = semipadd2pop_gt_data$Y2,
#' #'                                            Z2 = semipadd2pop_gt_data$Z2,
#' #'                                            Se2 = semipadd2pop_gt_data$Se2,
#' #'                                            Sp2 = semipadd2pop_gt_data$Sp2,
#' #'                                            X2 = semipadd2pop_gt_data$X2,
#' #'                                            nonparm2 = semipadd2pop_gt_data$nonparm2,
#' #'                                            rho1 = 2,
#' #'                                            rho2 = 1,
#' #'                                            w1 = 1,
#' #'                                            w2 = 1,
#' #'                                            w = 1,
#' #'                                            nCom = 4,
#' #'                                            d1 = semipadd2pop_gt_data$nonparm1 * 15,
#' #'                                            d2 = semipadd2pop_gt_data$nonparm2 * 10,
#' #'                                            xi = 1,
#' #'                                            n.lambda = 3,
#' #'                                            n.eta = 3,
#' #'                                            lambda.min.ratio =.001,
#' #'                                            lambda.beta = 1,
#' #'                                            lambda.f = 1,
#' #'                                            eta.beta = 1,
#' #'                                            eta.f = 1,
#' #'                                            tol = 1e-2,
#' #'                                            maxiter = 500,
#' #'                                            report.prog = TRUE)
#' #'
#' #' plot_semipadd2pop_gt_cv(semipadd2pop_gt_cv.out,
#' #'                        true.functions = list(f1 = semipadd2pop_gt_data$f1,
#' #'                                              f2 = semipadd2pop_gt_data$f2,
#' #'                                              X1 = semipadd2pop_gt_data$X1,
#' #'                                              X2 = semipadd2pop_gt_data$X2)
#' #' )
#' #' @export
#' semipadd2pop_gt_cv_adapt <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop_gt function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
#'   grouplasso2pop_gt_cv_adapt.out <- grouplasso2pop_gt_cv_adapt(Y1 = Y1,
#'                                                                Z1 = Z1,
#'                                                                Se1 = Se1,
#'                                                                Sp1 = Sp1,
#'                                                                X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                Y2 = Y2,
#'                                                                Z2 = Z2,
#'                                                                Se2 = Se2,
#'                                                                Sp2 = Sp2,
#'                                                                X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                rho1 = rho1,
#'                                                                rho2 = rho2,
#'                                                                n.lambda = n.lambda,
#'                                                                n.eta = n.eta,
#'                                                                lambda.min.ratio = lambda.min.ratio,
#'                                                                lambda.max.ratio = lambda.max.ratio,
#'                                                                n.folds = n.folds,
#'                                                                w1 = grouplasso2pop_inputs$w1,
#'                                                                w2 = grouplasso2pop_inputs$w2,
#'                                                                w = grouplasso2pop_inputs$w,
#'                                                                AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                Com = grouplasso2pop_inputs$Com,
#'                                                                E.approx = E.approx,
#'                                                                tol = tol,
#'                                                                maxiter = maxiter,
#'                                                                report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   f1.hat.folds <- vector("list",n.lambda)
#'   f2.hat.folds <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'     
#'     f1.hat.folds[[l]] <- vector("list",n.eta)
#'     f2.hat.folds[[l]] <- vector("list",n.eta)
#'     
#'     for( k in 1:n.eta)
#'     {
#'       
#'       f1.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       f2.hat.folds[[l]][[k]] <- vector("list",n.folds)
#'       
#'     }
#'     
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                nonparm1 = nonparm1,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                b1 = grouplasso2pop_gt_cv_adapt.out$b1.arr[,l,k],
#'                                                                X2 = X2,
#'                                                                nonparm2 = nonparm2,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                b2 = grouplasso2pop_gt_cv_adapt.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipadd2pop_gt_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipadd2pop_gt_fitted$beta2.hat
#'       
#'       # for( fold in 1:n.folds)
#'       # {
#'       # 
#'       #   semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'       #                                                           nonparm1 = nonparm1,
#'       #                                                           groups1 = grouplasso2pop_inputs$groups1,
#'       #                                                           knots.list1 = grouplasso2pop_inputs$knots.list1,
#'       #                                                           emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'       #                                                           QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'       #                                                           b1 = grouplasso2pop_gt_cv_adapt.out$b1.folds.arr[,l,k,fold],
#'       #                                                           X2 = X2,
#'       #                                                           nonparm2 = nonparm2,
#'       #                                                           groups2 = grouplasso2pop_inputs$groups2,
#'       #                                                           knots.list2 = grouplasso2pop_inputs$knots.list2,
#'       #                                                           emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'       #                                                           QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'       #                                                           b2 = grouplasso2pop_gt_cv_adapt.out$b2.folds.arr[,l,k,fold])
#'       # 
#'       #   f1.hat.folds[[l]][[k]][[fold]] <- semipadd2pop_gt_fitted$f1.hat
#'       #   f2.hat.folds[[l]][[k]][[fold]] <- semipadd2pop_gt_fitted$f2.hat
#'       # 
#'       # }
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.folds = f1.hat.folds,
#'                   f2.hat.folds = f2.hat.folds,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   P1.hat = grouplasso2pop_gt_cv_adapt.out$P1.arr,
#'                   P2.hat = grouplasso2pop_gt_cv_adapt.out$P2.arr,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = grouplasso2pop_gt_cv_adapt.out$lambda.seq,
#'                   eta.seq = grouplasso2pop_gt_cv_adapt.out$eta.seq,
#'                   lambda.initial.fit = grouplasso2pop_gt_cv_adapt.out$lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_gt_cv_adapt.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_gt_cv_adapt.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_gt_cv_adapt.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_gt_cv_adapt.out$w1,
#'                   w2 = grouplasso2pop_gt_cv_adapt.out$w2,
#'                   w = grouplasso2pop_gt_cv_adapt.out$w,
#'                   iterations = grouplasso2pop_gt_cv_adapt.out$iterations)
#'   
#'   class(output) <- "semipadd2pop_gt_cv"
#'   
#'   return(output)
#'   
#' }
#' 
#' 
#' 
#' #' Compute semiparametric regression model on group testing data with 2 data sets while penalizing dissimilarity
#' #'
#' #' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se1 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp1 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' #' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' #' @param Y2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Z2 Group testing output for data set 2 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' #' @param Se2 A vector of testing sensitivities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param Sp2 A vector of testing specificities, where the first element is the
#' #'      testing specificity for pools and the second entry is the
#' #'      test specificity for individual testing, if applicable.
#' #' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' #' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
#' #' @param rho1 weight placed on the first data set
#' #' @param rho2 weight placed on the second data set
#' #' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' #' @param w2 covariate-specific weights for different penalization among covariates in data set 2
#' #' @param w covariate-specific weights for different penalization toward similarity for different covariates
#' #' @param nCom the number of covariates to be treated as common between the two data sets: these must be arranged in the first \code{nCom} columns of the matrices \code{X1} and \code{X2} after the column of ones corresponding to the intercept.
#' #' @param d1 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param d2 vector giving the dimensions the B-spline bases to be used when fitting the nonparametric effects. If a scalar is given, this dimension is used for all nonparametric effects.
#' #' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' #' @param n.lambda the number of lambda values with which to make the grid
#' #' @param n.eta the number of eta values with which to make the grid
#' #' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' #' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' #' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' #' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' #' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' #' @param E.approx a logical indicating whether the conditional expectations in the E-step should be computed approximately or exactly.
#' #' @param tol a convergence criterion
#' #' @param maxiter the maximum allowed number of iterations (EM steps)
#' #' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' #' @return Returns the estimator of the semiparametric additive model with group testing data
#' #' @export
#' semipadd2pop_gt_cv_adapt_fixedgrid <- function(Y1,Z1,Se1,Sp1,X1,nonparm1,Y2,Z2,Se2,Sp2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.seq,eta.seq,lambda.initial.fit,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,E.approx = FALSE,tol=1e-3,maxiter = 1000,report.prog = FALSE)
#' {
#'   
#'   # prepare input for grouplasso2pop_gt function
#'   grouplasso2pop_inputs <- semipadd2pop_to_grouplasso2pop(X1 = X1,
#'                                                           nonparm1 = nonparm1,
#'                                                           X2 = X2,
#'                                                           nonparm2 = nonparm2,
#'                                                           nCom = nCom,
#'                                                           d1 = d1,
#'                                                           d2 = d2,
#'                                                           xi = xi,
#'                                                           w1 = w1,
#'                                                           w2 = w2,
#'                                                           w = w,
#'                                                           lambda.beta = lambda.beta,
#'                                                           lambda.f = lambda.f,
#'                                                           eta.beta = eta.beta,
#'                                                           eta.f = eta.f)
#'   
#'   # get group lasso estimators over a grid of lambda and eta values on the full data and CV folds
#'   grouplasso2pop_gt_cv_adapt_fixedgrid.out <- grouplasso2pop_gt_cv_adapt_fixedgrid(Y1 = Y1,
#'                                                                                    Z1 = Z1,
#'                                                                                    Se1 = Se1,
#'                                                                                    Sp1 = Sp1,
#'                                                                                    X1 = grouplasso2pop_inputs$DD1.tilde,
#'                                                                                    groups1 = grouplasso2pop_inputs$groups1,
#'                                                                                    Y2 = Y2,
#'                                                                                    Z2 = Z2,
#'                                                                                    Se2 = Se2,
#'                                                                                    Sp2 = Sp2,
#'                                                                                    X2 = grouplasso2pop_inputs$DD2.tilde,
#'                                                                                    groups2 = grouplasso2pop_inputs$groups2,
#'                                                                                    rho1 = rho1,
#'                                                                                    rho2 = rho2,
#'                                                                                    lambda.seq = lambda.seq,
#'                                                                                    eta.seq = eta.seq,
#'                                                                                    lambda.initial.fit = lambda.initial.fit,
#'                                                                                    n.folds = n.folds,
#'                                                                                    w1 = grouplasso2pop_inputs$w1,
#'                                                                                    w2 = grouplasso2pop_inputs$w2,
#'                                                                                    w = grouplasso2pop_inputs$w,
#'                                                                                    AA1 = grouplasso2pop_inputs$AA1.tilde,
#'                                                                                    AA2 = grouplasso2pop_inputs$AA2.tilde,
#'                                                                                    Com = grouplasso2pop_inputs$Com,
#'                                                                                    E.approx = E.approx,
#'                                                                                    tol = tol,
#'                                                                                    maxiter = maxiter,
#'                                                                                    report.prog = report.prog)
#'   
#'   # get matrices of the fitted functions evaluated at the design points
#'   
#'   n.lambda <- length(lambda.seq)
#'   n.eta <- length(eta.seq)
#'   
#'   f1.hat.design <- array(0,dim=c(nrow(X1),ncol(X1),n.lambda,n.eta))
#'   f2.hat.design <- array(0,dim=c(nrow(X2),ncol(X2),n.lambda,n.eta))
#'   
#'   beta1.hat <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
#'   beta2.hat <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
#'   
#'   f1.hat <- vector("list",n.lambda)
#'   f2.hat <- vector("list",n.lambda)
#'   
#'   for(l in 1:n.lambda)
#'   {
#'     
#'     f1.hat[[l]] <- vector("list",n.eta)
#'     f2.hat[[l]] <- vector("list",n.eta)
#'     
#'   }
#'   
#'   for(l in 1:n.lambda)
#'     for(k in 1:n.eta)
#'     {
#'       
#'       semipadd2pop_gt_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
#'                                                                nonparm1 = nonparm1,
#'                                                                groups1 = grouplasso2pop_inputs$groups1,
#'                                                                knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                                                                emp.cent1 = grouplasso2pop_inputs$emp.cent1,
#'                                                                QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
#'                                                                b1 = grouplasso2pop_gt_cv_adapt_fixedgrid.out$b1.arr[,l,k],
#'                                                                X2 = X2,
#'                                                                nonparm2 = nonparm2,
#'                                                                groups2 = grouplasso2pop_inputs$groups2,
#'                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
#'                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
#'                                                                b2 = grouplasso2pop_gt_cv_adapt_fixedgrid.out$b2.arr[,l,k])
#'       
#'       f1.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f1.hat
#'       f2.hat[[l]][[k]] <- semipadd2pop_gt_fitted$f2.hat
#'       
#'       f1.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f1.hat.design
#'       f2.hat.design[,,l,k] <- semipadd2pop_gt_fitted$f2.hat.design
#'       
#'       beta1.hat[,l,k] <- semipadd2pop_gt_fitted$beta1.hat
#'       beta2.hat[,l,k] <- semipadd2pop_gt_fitted$beta2.hat
#'       
#'       
#'     }
#'   
#'   # prepare output
#'   output <- list( f1.hat = f1.hat,
#'                   f2.hat = f2.hat,
#'                   f1.hat.design = f1.hat.design,
#'                   f2.hat.design = f2.hat.design,
#'                   P1.hat = grouplasso2pop_gt_cv_adapt_fixedgrid.out$P1.arr,
#'                   P2.hat = grouplasso2pop_gt_cv_adapt_fixedgrid.out$P2.arr,
#'                   beta1.hat = beta1.hat,
#'                   beta2.hat = beta2.hat,
#'                   rho1 = rho1,
#'                   rho2 = rho2,
#'                   nonparm1 = nonparm1,
#'                   nonparm2 = nonparm2,
#'                   d1 = d1,
#'                   d2 = d2,
#'                   xi = xi,
#'                   knots.list1 = grouplasso2pop_inputs$knots.list1,
#'                   knots.list2 = grouplasso2pop_inputs$knots.list2,
#'                   lambda.beta = lambda.beta,
#'                   lambda.f = lambda.f,
#'                   eta.beta = eta.beta,
#'                   eta.f = eta.f,
#'                   Com = grouplasso2pop_inputs$Com,
#'                   n.lambda = n.lambda,
#'                   n.eta = n.eta,
#'                   n.folds = n.folds,
#'                   lambda.seq = lambda.seq,
#'                   eta.seq = eta.seq,
#'                   lambda.initial.fit = lambda.initial.fit,
#'                   which.lambda.cv = grouplasso2pop_gt_cv_adapt_fixedgrid.out$which.lambda.cv,
#'                   which.eta.cv = grouplasso2pop_gt_cv_adapt_fixedgrid.out$which.eta.cv,
#'                   which.lambda.cv.under.zero.eta = grouplasso2pop_gt_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
#'                   w1 = grouplasso2pop_gt_cv_adapt_fixedgrid.out$w1,
#'                   w2 = grouplasso2pop_gt_cv_adapt_fixedgrid.out$w2,
#'                   w = grouplasso2pop_gt_cv_adapt_fixedgrid.out$w,
#'                   iterations = grouplasso2pop_gt_cv_adapt_fixedgrid.out$iterations)
#'   
#'   class(output) <- "semipadd2pop_gt_cv"
#'   
#'   return(output)
#'   
#' }
#' 
