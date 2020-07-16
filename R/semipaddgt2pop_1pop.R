#' Compute the objective function of the 1-population group lasso problem with a binary response
#'
#' @param beta1 the vector of coefficients for data set 1
#' @param Y1 the response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param lambda the level of sparsity penalization
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @return Returns the value of the objective function for the 2-population group lasso problem
#' @export
grouplasso1pop_logreg_obj <- function(beta1,Y1,X1,groups1,lambda,w1)
{

  q1 <- length(unique(groups1))
  n1 <- nrow(X1)

  P1 <- logit(X1 %*% beta1)

  neg2LL1 <- - 2 * sum(  Y1 * log(P1) + (1 - Y1)*log(1-P1))

  beta1.wl2l1 <- 0
  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    beta1.wl2l1  <- beta1.wl2l1 + w1[j] * sqrt(sum( beta1[ind]^2 ))
  }

  pen1 <- lambda * beta1.wl2l1

  val <- neg2LL1 + pen1

  return(val)

}
#' Minimize the objective function of the 1-population group lasso problem with a binary response
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
#' @param lambda the level of sparsity penalization
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the minimizer of the 1-population group lasso objective function for the two data sets.
#'
#' @examples
#' # # generate data
#' n1 <- 400
#' d1 <- c(1,1,3,4)
#' q1 <- length(d1)
#' X1 <- matrix(rnorm(n1*sum(d1)),n1,sum(d1))
#' groups1 <- numeric() ; for(j in 1:q1){ groups1 <- c(groups1,rep(j,d1[j])) }
#' beta1 <- rnorm(ncol(X1))
#' Y1 <- rbinom(n1,1,logit(X1 %*% beta1))
#'
#' # set tuning parameters
#' lambda <- 10
#' w1 <- rexp(q1,2)
#'
#' # fit grouplasso1pop_logreg estimator
#' grouplasso1pop_logreg.out <- grouplasso1pop_logreg(Y1,X1,groups1,lambda,w1,tol=1e-4,max.iter=500,return_obj=TRUE)
#'
#' beta1.hat <- grouplasso1pop_logreg.out$beta1.hat
#' obj.val <- grouplasso1pop_logreg.out$obj.val
#'
#' # look at results
#' beta1.hat
#' plot(obj.val)
#' any(diff(obj.val) > 0) == FALSE # does the objective function decrease after every iteration?  YESSS!!!!!!
#' @export
grouplasso1pop_logreg <- function(Y1,X1,groups1,lambda,w1,tol=1e-4,max.iter=500,return_obj=FALSE,init = NULL)
{

  # initialize estimators and convergence criteria
  q1 <- ncol(X1)
  n1 <- nrow(X1)

  if( length(init) == 0){

    beta1.hat1 <- rep(0,q1)

  } else {

    beta1.hat1 <- init

  }


  conv <- 1
  iter <- 0

  got.eigen1 <- numeric(q1)
  eigen1 <- vector("list", q1)


  d1 <- table(groups1)
  singleton.grps1 <- which(d1 == 1)
  nonsingleton.grps1 <- which(d1 != 1)
  obj.val <- numeric()

  while( conv > tol & iter < max.iter)
  {

    beta1.hat0 <- beta1.hat1


    # go through the groups of size 1 of the first data set
    for(j in singleton.grps1)
    {

      P1 <- logit( X1 %*% beta1.hat1)

      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind1 <- which(groups1 == j)
      X1j <- X1[,ind1,drop=FALSE]
      X1j.tilde <- diag(sqrt(P1*(1-P1))) %*% X1j
      Z1j.tilde <- diag(sqrt(P1*(1-P1))) %*% ( X1j %*% beta1.hat1[ind1] + diag(1/(P1*(1-P1))) %*% ( Y1 - P1 ) )
      sx1j <- sum(X1j.tilde^2)

      sx1z1j <- sum(X1j.tilde * Z1j.tilde)
      beta1.hat1[ind1] <- SoftThresh(sx1z1j,lambda*w1[j]/2)/sx1j

    }


    # go through the groups of size more than 1 of the first data set
    for(j in nonsingleton.grps1)
    {

      P1 <- logit( X1 %*% beta1.hat1)

      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind1 <- which(groups1 == j)
      X1j <- X1[,ind1,drop=FALSE]

      # fixed Hessian approach for the groups
      X1j.tilde <- (1/2)*X1j
      Z1j.tilde <- X1j.tilde %*% beta1.hat1[ind1] + 2 * ( Y1 - P1 )

      x1z1j.norm <- sqrt(sum(  (t(X1j.tilde) %*% Z1j.tilde)^2 ))

      if( x1z1j.norm < lambda * w1[j]/2)
      {

        beta1.hat1[ind1] <- 0

      } else {

        # compute FD update
        if(got.eigen1[j]==0){

          LtL <- t(X1j.tilde) %*% X1j.tilde
          eigen1[[j]] <- eigen(LtL)
          got.eigen1[j] <- 1

        }

        h <- Z1j.tilde
        L <- X1j.tilde

        beta1.hat1[ind1] <- FoygelDrton(h,L,lambda*w1[j]/2,eigen1[[j]]$values,t(eigen1[[j]]$vectors))

      }

    }

    conv <- max(abs(c(beta1.hat1 - beta1.hat0)))
    iter <- iter + 1

    if(return_obj == TRUE)
    {
      obj.val[iter] <- grouplasso1pop_logreg_obj(beta1.hat1,Y1,X1,groups1,lambda,w1)
    }

  }

  beta1.hat <- beta1.hat1

  output <- list(beta1.hat = beta1.hat,
                 obj.val = obj.val,
                 iter = iter)

  return(output)

}

#' Compute semiparametric binary-response regression model with 1 data set
#'
#' @param Y1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z1 Group testing output for data set 1 in the format as output by one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se1 A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param Sp1 A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the
#'      test specificity for individual testing, if applicable.
#' @param XX1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param w1 covariate-specific weights for different penalization among covariates in data set 1
#' @param d the dimension of the B-spline basis to be used when fitting the nonparametric effects
#' @param xi a tuning parameter governing the smoothness of the nonparametric estimates
#' @param lambda.beta the level of sparsity penalization for the parametric effects
#' @param lambda.f the level of sparsity penalization for the nonparametric effects
#' @param tol a convergence criterion
#' @param max.iter the maximum allowed number of iterations (EM steps)
#' @param report.prog a logical. If \code{TRUE} then the number of inner loops required to complete the M step of the EM algorithm are returned after each EM step.
#' @return Returns the estimator of the semiparametric additive model with group testing data
#'
#' @examples
#' # generate data set
#' n1 <- 200
#' p1 <- 6
#' q1 <- 4
#' zeta1 <- 3/20
#' zeta2 <- 10/20
#' W1 <- cbind(corrBern(n1,probs=c(1:p1)/(2*p1),Rho = zeta1^abs(outer(1:p1,1:p1,"-"))))
#' X1 <- (corrUnif(n1,Rho = zeta2^abs(outer(1:q1,1:q1,"-")))-.5)*5
#'
#' XX1 <- cbind(1,W1[,c(1,2)],X1[,c(1,2)],W1[,-c(1,2)],X1[,-c(1,2)])
#' nonparm1 <- c(0,rep(0,2),rep(1,2),rep(0,p1 - 2),rep(1,q1 - 2))
#' pp1 <- ncol(XX1)
#'
#' # set up the true functions
#' f1 <- vector("list",pp1)
#' f1[[1]] <- function(x){0} # intercept
#' f1[[2]] <- function(x){x*2}
#' f1[[3]] <- function(x){x*0}
#' f1[[4]] <- function(x){-sin(x*2) }
#' f1[[5]] <- function(x){x}
#' f1[[6]] <- function(x){x*2}
#' f1[[7]] <- function(x){x*0}
#' f1[[8]] <- function(x){0*x}
#' f1[[9]] <- function(x){0*x}
#' f1[[10]] <- function(x){0*x}
#' f1[[11]] <- function(x){exp(-x)-2/5*sinh(5/2)}
#'
#' f1.design <- matrix(0,n1,pp1)
#' for(j in 1:pp1){
#'   if(nonparm1[j]==1){
#'     f1.design[,j] <- f1[[j]](XX1[,j]) - mean(f1[[j]](XX1[,j]))
#'   } else {
#'     f1.design[,j] <- f1[[j]](XX1[,j])
#'   }
#' }
#'
#' Y1.true <- rbinom(n1,1,logit(apply(f1.design,1,sum)))
#'
#' Se1 <- c(.98,.96)
#' Sp1 <- c(.97,.99)
#' assay1.out <- dorfman.assay.gen(Y1.true,Se1,Sp1,cj=4)
#' Z1 <- assay1.out$Z
#' Y1 <- assay1.out$Y
#'
#' semipaddgt1pop.out <- semipaddgt1pop_logreg(Y1,
#'                                      Z1,
#'                                      Se1,
#'                                      Sp1,
#'                                      XX1,
#'                                      nonparm1,
#'                                      w1=1,
#'                                      d=20,
#'                                      xi=.5,
#'                                      lambda.beta=2,
#'                                      lambda.f=1,
#'                                      tol=1e-3,
#'                                      max.iter = 50,
#'                                      report.prog = TRUE)
#'
#' f1.hat.design <- semipaddgt1pop.out$f1.hat.design
#' f1.hat <- semipaddgt1pop.out$f1.hat
#' knots.list1 <- semipaddgt1pop.out$knots.list1
#' nonparm1 <- semipaddgt1pop.out$nonparm1
#'
#' # number of iterations
#' semipaddgt1pop.out$inner.iter
#'
#' # plot results
#' par(mfrow=c(3,4),mar=c(0,0,0,0))
#' for( j in 2:pp1){
#'
#'   plot(NA,ylim = range(f1.hat.design),xlim=range(XX1[,j]),xaxt="n",yaxt="n")
#'   if(nonparm1[j]==1) abline(v=knots.list1[[j]],col=rgb(0,0,0,.15))
#'
#'   plot(f1.hat[[j]],min(XX1[,j]),max(XX1[,j]),add=TRUE)
#'   x1.seq <- seq(min(XX1[,j]),max(XX1[,j]),length=300)
#'   f1.cent.seq <- f1[[j]](x1.seq) - mean(f1[[j]](XX1[,j]))
#'   lines(f1.cent.seq ~ x1.seq,lty=2)
#' }
#' @export
semipaddgt1pop_logreg <- function(Y1,Z1,Se1,Sp1,XX1,nonparm1,w1,d,xi,lambda.beta,lambda.f,tol=1e-3,max.iter=500,report.prog=FALSE)
{

  ww1 <- ((1-nonparm1) + nonparm1 * lambda.f/lambda.beta) * w1
  ww1[1] <- 0 # do not penalize intercept

  # prepare input for grouplasso2pop function
  n1 <- nrow(XX1)
  pp1 <- ncol(XX1)

  ##### For first data set
  DD1.tilde <- matrix(NA,n1,0)
  groups1 <- numeric()
  QQ1.inv <- vector("list",length=pp1)
  knots.list1 <- vector("list",length=pp1)
  emp.cent1 <- vector("list",length=pp1)

  for( j in 1:pp1 )
  {

    if(nonparm1[j] == 0){

      DD1.tilde <- cbind(DD1.tilde,XX1[,j])
      groups1 <- c(groups1,j)

    } else {


      int.knots <- quantile(XX1[,j],seq(0,1,length=d-2+1)) # add one, so that one can be removed after centering to restore full-rank.
      boundary.knots <- range(int.knots)
      all.knots <- sort(c(rep(boundary.knots,3),int.knots))
      knots.list1[[j]] <- all.knots

      B1j <- spline.des(all.knots,XX1[,j],ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
      emp.cent1[[j]] <- apply(B1j,2,mean)
      B1j.cent <- B1j - matrix(emp.cent1[[j]],n1,d,byrow=TRUE)

      # construct matrix in which l2 norm of function is a quadratic form
      M <- t(B1j.cent) %*% B1j.cent / n1

      # construct matrix in which 2nd derivative penalty is a quadratic form
      R <- matrix(NA,d+1,d+1)
      dsq_bspline.mat <- spline.des(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)$design
      for(k in 1:(d+1))
        for(l in 1:(d+1))
        {

          pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
          h <- diff(int.knots)
          R[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.

        }

      Q <-  chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
      Q.inv <- solve(Q)

      QQ1.inv[[j]] <- Q.inv

      # construct transformed basis function matrices
      DD1.tilde <- cbind(DD1.tilde, B1j.cent %*% Q.inv)
      groups1 <- c(groups1,rep(j,d))

    }

  }

  Y1.diag <- pull.diagnoses(Z1,Y1)

  if( length(Y1.diag) != 0){

    grouplasso1pop_logreg.out <- grouplasso1pop_logreg(Y1.diag,
                                                       X1 = DD1.tilde,
                                                       groups1,
                                                       lambda = lambda.beta,
                                                       w1 = ww1,
                                                       tol=tol,
                                                       max.iter=max.iter,
                                                       return_obj=FALSE)

    beta1.hat.diag <- grouplasso1pop_logreg.out$beta1.hat # should I return the fitted functions according to the diagnoses?
    beta1.hat1 <- beta1.hat.diag


  } else {

    beta1.hat1 <- rep(0,ncol(DD1.tilde))
  }

  ###### Do the EM-algorithm with penalized updates
  conv <- 1
  iter <- 0
  inner.iter <- numeric()
  while( conv > tol & iter < max.iter)
  {

    beta1.hat0 <- beta1.hat1

    # E-step: compute the conditional expectations for the true disease statuses
    EY1 <- EY.exact(Z1,Y1,X=DD1.tilde,b=beta1.hat1,Se1,Sp1)

    # update initial values
    init <- beta1.hat1

    # M-step: maximize the objective function with conditional expectations substituted
    grouplasso1pop_logreg.out <- grouplasso1pop_logreg(EY1,
                                                       X1 = DD1.tilde,
                                                       groups1,
                                                       lambda = lambda.beta,
                                                       w1 = ww1,
                                                       tol=tol,
                                                       max.iter = max.iter,
                                                       return_obj=FALSE,
                                                       init = init )

    beta1.hat1 <- grouplasso1pop_logreg.out$beta1.hat

    conv <- max(abs(c(beta1.hat1 - beta1.hat0)))
    iter <- iter + 1
    inner.iter[iter] <- grouplasso1pop_logreg.out$iter
    if(report.prog) print(grouplasso1pop_logreg.out$iter)

  }

  b1 <- beta1.hat1


  # store fitted functions on data set 1 in a list
  f1.hat <- vector("list",pp1)
  f1.hat.design <- matrix(0,n1,pp1)

  f1.hat[[1]] <- eval( parse( text= paste("function(x){",paste(b1[1])," }")))
  f1.hat.design[,1] <- b1[1]
  for(j in 2:pp1)
  {

    if(nonparm1[j] == 0)
    {

      ind <- which(groups1 == j)
      f1.hat[[j]] <- eval( parse( text= paste("function(x){ x * ",paste(b1[ind])," }")))

    } else {

      ind <- which(groups1 == j)

      Gamma1.hat <- QQ1.inv[[j]] %*% b1[ind]
      f1.hat[[j]] <- eval(parse(text=paste("function(x)
                                           {

                                           x <- round(x,10)
                                           x.mat <- spline.des(",paste("c(",paste(round(knots.list1[[j]],10),collapse=","),")",sep=""),",x,ord=4,derivs=0,outer.ok=TRUE)$design[,-1]
                                           x.mat.cent <- x.mat - matrix(",paste("c(",paste(emp.cent1[[j]],collapse=","),"),length(x),",d,sep=""),",byrow=TRUE)
                                           f.hat <- as.numeric(x.mat.cent %*% ",paste("c(",paste(Gamma1.hat,collapse=","),")",sep=""),")

                                           return(f.hat)

                                           }"
        			)))

    }

    f1.hat.design[,j] <- f1.hat[[j]](XX1[,j])

}



  # collect output
  output <- list( f1.hat = f1.hat,
                  f1.hat.design = f1.hat.design,
                  nonparm1 = nonparm1,
                  d = d,
                  xi = xi,
                  knots.list1 = knots.list1,
                  lambda.beta = lambda.beta,
                  lambda.f = lambda.f,
                  inner.iter = inner.iter
  )

  return(output)

}
