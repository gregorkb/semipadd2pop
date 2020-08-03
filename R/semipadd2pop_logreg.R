#' Compute the objective function of the 2-population group lasso problem with a binary response
#'
#' @param beta1 the vector of coefficients for data set 1
#' @param beta2 the vector of coefficients for data set 2
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
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
#' @return Returns the value of the objective function for the 2-population group lasso problem
#' @export
grouplasso2pop_logreg_obj <- function(beta1,beta2,Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda,eta,w1,w2,w,AA1,AA2,Com)
{

  q1 <- length(unique(groups1))
  q2 <- length(unique(groups2))

  n1 <- nrow(X1)
  n2 <- nrow(X2)

  P1 <- logit(X1 %*% beta1)
  P2 <- logit(X2 %*% beta2)

  # add weights to objective function:
  neg2LL1 <- - 2 * rho1 / n1 * sum(  Y1 * log(P1) + (1 - Y1)*log(1-P1) )
  neg2LL2 <- - 2 * rho2 / n2 * sum(  Y2 * log(P2) + (1 - Y2)*log(1-P2) )

  beta1.wl2l1 <- 0
  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    beta1.wl2l1  <- beta1.wl2l1 + w1[j] * sqrt(sum( beta1[ind]^2 ))
  }

  beta2.wl2l1 <- 0
  for(j in 1:q2)
  {
    ind <- which(groups2 == j)
    beta2.wl2l1  <- beta2.wl2l1 + w2[j] * sqrt(sum( beta2[ind]^2 ))
  }

  beta1beta2.wl2 <- 0
  for(j in Com)
  {
    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)
    beta1beta2.wl2  <- beta1beta2.wl2 + w[j] * mean( (AA1[[j]] %*% beta1[ind1] - AA2[[j]] %*% beta2[ind2] )^2 )
  }

  pen1 <- lambda * (beta1.wl2l1 + (rho2/rho1) * beta2.wl2l1 )

  pen2 <- eta * beta1beta2.wl2

  val <- neg2LL1 + neg2LL2 + pen1 + pen2

  return(val)

}


#' Minimize the objective function of the 2-population group lasso problem with a binary response
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector of integers indicating to which group each covariate in data set 1 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda the level of sparsity penalization
#' @param eta the level of penalization towards model similarity
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param eigen1 a list of eigen info on groups from data set 1
#' @param eigen2 a list of eigen info on groups from data set 2
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @param init a list with vectors \code{beta1} and \code{beta2} to serve as initial values for the parameters
#' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
#'
#' @examples
#' # generate data
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 400,
#'                                                              n2 = 600)
#'
#' # fit grouplasso2pop estimator
#' grouplasso2pop_logreg_R.out <- grouplasso2pop_logreg_R(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                    X1 = grouplasso2pop_logreg_data$X1,
#'                                                    groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                    Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                    X2 = grouplasso2pop_logreg_data$X2,
#'                                                    groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                    rho1 = 2,
#'                                                    rho2 = 1,
#'                                                    lambda = 1,
#'                                                    eta = 1,
#'                                                    w1 = grouplasso2pop_logreg_data$w1,
#'                                                    w2 = grouplasso2pop_logreg_data$w2,
#'                                                    w = grouplasso2pop_logreg_data$w,
#'                                                    AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                    AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                    Com = grouplasso2pop_logreg_data$Com,
#'                                                    tol = 1e-4,
#'                                                    maxiter = 500,
#'                                                    plot_obj = TRUE)
#' @export
grouplasso2pop_logreg_R <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda,eta,w1,w2,w,AA1,AA2,Com,tol=1e-4,maxiter=500,plot_obj=FALSE,init = NULL)
{
  
  # initialize estimators and convergence criteria
  q1 <- ncol(X1)
  q2 <- ncol(X2)
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  if( length(init) == 0){
    
    beta1.hat1 <- rep(0,q1)
    beta2.hat1 <- rep(0,q2)
    
  } else {
    
    beta1.hat1 <- init$beta1
    beta2.hat1 <- init$beta2
    
  }
  
  conv <- 1
  iter <- 0
  
  d1 <- table(groups1)
  d2 <- table(groups2)
  singleton.grps1 <- which(d1 == 1)
  nonsingleton.grps1 <- which(d1 != 1)
  singleton.grps2 <- which(d2 == 1)
  nonsingleton.grps2 <- which(d2 != 1)
  
  got.eigen1 <- numeric(length(d1))
  got.eigen2 <- numeric(length(d2))
  eigen1 <- vector("list", length(d1))
  eigen2 <- vector("list", length(d2))
  
  obj.val <- numeric()
  
  lambda1 <- lambda * n1 / rho1
  eta1 <- eta * n1 / rho1
  
  lambda2 <- lambda * n2 / rho1 # because it should be rho2/rho1 * lambda * n2 / rho2.
  eta2 <- eta * n2 / rho2
  
  while( conv > tol & iter < maxiter)
  {
    
    beta1.hat0 <- beta1.hat1
    beta2.hat0 <- beta2.hat1
    
    # go through the groups of size 1 of the first data set
    
    for(j in singleton.grps1)
    {
      
      P1 <- logit( X1 %*% beta1.hat1)
      if( any(P1  > 1 - 1e-10) | any(P1 < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # fixed Hessian
      ind1 <- which(groups1 == j)
      X1j <- X1[,ind1,drop=FALSE]
      X1j.tilde <- (1/2) * X1j
      Z1j.tilde <- (1/2) * ( X1j %*% beta1.hat1[ind1] + 4 * ( Y1 - P1 ) )
      sx1j <- sum(X1j.tilde^2)
      
      if(j %in% Com){
        
        ind2 <- which(groups2 == j)
        sx1z1j.AAb2 <- sum(X1j.tilde * Z1j.tilde) + eta1*w[j]*t(AA1[[j]]) %*% AA2[[j]] %*% beta2.hat1[ind2]
        beta1.hat1[ind1] <- SoftThresh_R(sx1z1j.AAb2,lambda1*w1[j]/2)/(sx1j + eta1*w[j]*t(AA1[[j]]) %*% AA1[[j]])
        
        if(is.na(sx1z1j.AAb2)) print(j)
        
      } else {
        
        sx1z1j <- sum(X1j.tilde * Z1j.tilde)
        beta1.hat1[ind1] <- SoftThresh_R(sx1z1j,lambda1*w1[j]/2)/sx1j
        
      }
      
    }
    
    # go through the groups of size more than 1 of the first data set
    for(j in nonsingleton.grps1)
    {
      
      P1 <- logit( X1 %*% beta1.hat1)
      if( any(P1  > 1 - 1e-10) | any(P1 < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind1 <- which(groups1 == j)
      X1j <- X1[,ind1,drop=FALSE]
      
      # fixed Hessian approach for the groups
      X1j.tilde <- (1/2)*X1j
      Z1j.tilde <- (1/2)*(X1j %*% beta1.hat1[ind1] + 4*(Y1 - P1))
      
      # if j in Com
      if(j %in% Com)
      {
        
        ind2 <- which(groups2 == j)
        x1z1j.AAb2 <- t(X1j.tilde) %*% Z1j.tilde + eta1*w[j]*t(AA1[[j]]) %*% AA2[[j]] %*% beta2.hat1[ind2]
        x1z1j.AAb2.norm <- sqrt(sum( x1z1j.AAb2^2 ))
        
        if(x1z1j.AAb2.norm < lambda1 * w1[j]/2){
          
          beta1.hat1[ind1] <- 0
          
        } else {
          
          # compute modified DF update
          if(got.eigen1[j]==0){
            
            LtL <- t(X1j.tilde) %*% X1j.tilde + eta1 * w[j] * t(AA1[[j]])%*% AA1[[j]]
            eigen1[[j]] <- eigen( LtL )
            eigen1[[j]][["LtLchol"]] <- chol(LtL)
            got.eigen1[j] <- 1
            
          }
          
          L <- eigen1[[j]]$LtLchol
          h <- solve(t(L)) %*% x1z1j.AAb2
          
          beta1.hat1[ind1] <- FoygelDrton_R(h,L,lambda1*w1[j]/2,eigen1[[j]]$values,t(eigen1[[j]]$vectors))
          
        }
        
        # if j not in Com
      } else {
        
        x1z1j.norm <- sqrt(sum(  (t(X1j.tilde) %*% Z1j.tilde)^2 ))
        
        if( x1z1j.norm < lambda1 * w1[j]/2)
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
          
          beta1.hat1[ind1] <- FoygelDrton_R(h,L,lambda1*w1[j]/2,eigen1[[j]]$values,t(eigen1[[j]]$vectors))
          
        }
        
      } # end if j not in Com
      
    } # end going through groups of first data set
    

    # go through the groups of size 1 of the second data set
    for(j in singleton.grps2)
    {
      
      P2 <- logit( X2 %*% beta2.hat1)
      if( any(P2  > 1 - 1e-10) | any(P2 < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind2 <- which(groups2 == j)
      X2j <- X2[,ind2,drop=FALSE]
      X2j.tilde <- (1/2)* X2j
      Z2j.tilde <- (1/2)* ( X2j %*% beta2.hat1[ind2] + 4 * ( Y2 - P2 ) )
      sx2j <- sum(X2j.tilde^2)
      
      if(j %in% Com){
        
        ind1 <- which(groups1 == j)
        sx2z2j.AAb2 <- sum(X2j.tilde * Z2j.tilde) + eta2*w[j]*t(AA2[[j]]) %*% AA1[[j]] %*% beta1.hat1[ind1]
        beta2.hat1[ind2] <- SoftThresh_R(sx2z2j.AAb2,lambda2*w2[j]/2)/(sx2j + eta2*w[j]*t(AA2[[j]]) %*% AA2[[j]])
        
      } else {
        
        sx2z2j <- sum(X2j.tilde * Z2j.tilde)
        beta2.hat1[ind2] <- SoftThresh_R(sx2z2j ,lambda2*w2[j]/2)/sx2j
        
      }
      
    }
    
    # go through the groups of size more than 1 of the second data set
    for(j in nonsingleton.grps2)
    {
      
      P2 <- logit( X2 %*% beta2.hat1)
      if( any(P2  > 1 - 1e-10) | any(P2 < 1e-10) ) stop("fitted probabilities diverging to 1 or 0")
      
      # set up iteratively re-weighted least-squares design matrix and response vector for coordinate j
      ind2 <- which(groups2 == j)
      X2j <- X2[,ind2,drop=FALSE]
      
      # fixed Hessian approach for the groups
      X2j.tilde <- (1/2) * X2j
      Z2j.tilde <- (1/2)*(X2j %*% beta2.hat1[ind2] + 4*(Y2 - P2))
      
      # if j in Com
      if(j %in% Com)
      {
        
        ind1 <- which(groups1 == j)
        x2z2j.AAb2 <- t(X2j.tilde) %*% Z2j.tilde + eta2*w[j]*t(AA2[[j]]) %*% AA1[[j]] %*% beta1.hat1[ind1]
        x2z2j.AAb2.norm <- sqrt(sum( x2z2j.AAb2^2 ))
        
        if(x2z2j.AAb2.norm < lambda2 * w2[j]/2)
        {
          
          beta2.hat1[ind2] <- 0
          
        } else {
          
          # compute modified DF update
          if(got.eigen2[j]==0){
            
            LtL <- t(X2j.tilde) %*% X2j.tilde + eta2 * w[j] * t(AA2[[j]])%*% AA2[[j]]
            eigen2[[j]] <- eigen( LtL )
            eigen2[[j]][["LtLchol"]] <- chol(LtL)
            got.eigen2[j] <- 1
            
          }
          
          L <- eigen2[[j]]$LtLchol
          h <- solve(t(L)) %*% x2z2j.AAb2
          
          beta2.hat1[ind2] <- FoygelDrton_R(h,L,lambda2*w2[j]/2,eigen2[[j]]$values,t(eigen2[[j]]$vectors))
          
        }
        
        # if j not in Com
      } else {
        
        x2z2j.norm <- sqrt(sum(  (t(X2j.tilde) %*% Z2j.tilde)^2 ))
        
        if( x2z2j.norm < lambda2 * w2[j]/2)
        {
          
          beta2.hat1[ind2] <- 0
          
        } else {
          
          # compute FD update
          if(got.eigen2[j]==0){
            
            LtL <- t(X2j.tilde) %*% X2j.tilde
            eigen2[[j]] <- eigen(LtL)
            got.eigen2[j] <- 1
            
          }
          
          h <- Z2j.tilde
          L <- X2j.tilde
          
          beta2.hat1[ind2] <- FoygelDrton_R(h,L,lambda2*w2[j]/2,eigen2[[j]]$values,t(eigen2[[j]]$vectors))
          
        }
        
      } # end if j not in Com
      
    } # end going through groups of second data set
    
    conv <- max(abs(c(beta1.hat1 - beta1.hat0, beta2.hat1 - beta2.hat0)))
    iter <- iter + 1
    
    if(plot_obj == TRUE)
    {
      obj.val[iter] <- grouplasso2pop_logreg_obj(beta1.hat1,beta2.hat1,Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda,eta,w1,w2,w,AA1,AA2,Com)
    }
    
  }
  
  if(plot_obj == TRUE)
  {
    plot(obj.val)
  }
  
  beta1.hat <- beta1.hat1
  beta2.hat <- beta2.hat1
  
  output <- list(beta1.hat = beta1.hat,
                 beta2.hat = beta2.hat,
                 obj.val = obj.val,
                 iter = iter)
  
  return(output)
  
}


#' Fit grouplasso2pop logistic regression estimator over a grid of lambda and eta values
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values desired
#' @param n.eta the number of eta values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#' @examples 
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 400,
#'                                                              n2 = 600)
#' 
#' grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                              X1 = grouplasso2pop_logreg_data$X1,
#'                                                              groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                              Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                              X2 = grouplasso2pop_logreg_data$X2,
#'                                                              groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                              rho1 = 2,
#'                                                              rho2 = 1,
#'                                                              n.lambda = 5,
#'                                                              n.eta = 5,
#'                                                              lambda.min.ratio = 0.01,
#'                                                              w1 = grouplasso2pop_logreg_data$w1,
#'                                                              w2 = grouplasso2pop_logreg_data$w2,
#'                                                              w = grouplasso2pop_logreg_data$w,
#'                                                              AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                              AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                              Com = grouplasso2pop_logreg_data$Com,
#'                                                              tol = 1e-3,
#'                                                              maxiter = 500,
#'                                                              report.prog = TRUE)
#' @export
grouplasso2pop_logreg_grid <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,lambda.max.ratio = 1,w1,w2,w,AA1,AA2,Com,tol=1e-4,maxiter=500,report.prog=FALSE)
{

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

    if(j <= q1)
    {
      norms1[j] <- sqrt(sum((t(X1[,ind1]) %*% (Y1 - mean(Y1)))^2)) / ( w1[j] * n1 / rho1)
    }
    if(j <= q2)
    {
      norms2[j] <- sqrt(sum((t(X2[,ind2]) %*% (Y2 - mean(Y2)))^2)) / (w2[j] * n2 / rho1)
    }

  }

  lambda.max <- 2 * max(norms1,norms2) # yes this is correct! This is the smallest value of lambda which sets all the non-intercept entries of beta1 and beta2 equal to zero.

  # make a lambda sequence
  largest.lambda <- lambda.max.ratio * lambda.max
  smallest.lambda <- lambda.min.ratio * lambda.max
  lambda.seq <- sort(c(exp(log(smallest.lambda) + ((n.lambda+1):1)/(n.lambda+1) * ((log(largest.lambda) - log(smallest.lambda)))))[-1])

  if(n.lambda == 1) lambda.seq <- smallest.lambda

  # make the eta sequence after fitting the model for the smallest value of lambda
  eta.seq <- numeric(n.eta)

  b1.arr <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  b2.arr <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
  
  P1.arr <- array(0,dim=c(nrow(X1),n.lambda,n.eta))
  P2.arr <- array(0,dim=c(nrow(X2),n.lambda,n.eta))

  iterations <- matrix(0,n.lambda*n.eta,3)
  colnames(iterations) <- c("lambda","eta","iter")
  step <- 0
  init <- list( beta1 = rep(0,q1),
                beta2 = rep(0,q2))
  for(l in 1:n.lambda){
    for(k in 1:n.eta){

      grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
                                                         rX1 = X1,
                                                         groups1 = groups1,
                                                         rY2 = Y2,
                                                         rX2 = X2,
                                                         groups2 = groups2,
                                                         rho1 = rho1,
                                                         rho2 = rho2,
                                                         lambda = lambda.seq[l],
                                                         eta = eta.seq[k],
                                                         w1 = w1,
                                                         w2 = w2,
                                                         w = w,
                                                         rAA1 = AA1,
                                                         rAA2 = AA2,
                                                         rCom = Com,
                                                         tol = tol,
                                                         maxiter = maxiter,
                                                         beta1_init = init$beta1,
                                                         beta2_init = init$beta2)

      b1 <- grouplasso2pop_logreg.out$beta1.hat
      b2 <- grouplasso2pop_logreg.out$beta2.hat

      if(l == 1 & k == 1){# define the eta sequence

        P1 <- logit(X1 %*% b1)
        P2 <- logit(X2 %*% b2)

        neg2LL1 <- - 2 * rho1 / n1 * sum(Y1*log(P1) + (1 - Y1)*log(1-P1))
        neg2LL2 <- - 2 * rho2 / n2 * sum(Y2*log(P2) + (1 - Y2)*log(1-P2))

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

      P1.arr[,l,k] <- logit(X1 %*% b1)
      P2.arr[,l,k] <- logit(X2 %*% b2)
      
      step <- step + 1
      iterations[step,] <- c(lambda.seq[l],eta.seq[k],grouplasso2pop_logreg.out$iter)

      if(report.prog == TRUE){

        print(c(l,k,grouplasso2pop_logreg.out$iter))

      }

    }

  }

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  P1.arr = P1.arr,
                  P2.arr = P2.arr,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = iterations)

  return(output)

}

#' Fit grouplasso2pop logistic regression estimator over a user-specified grid of lambda and eta values
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
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
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#' @export
grouplasso2pop_logreg_fixedgrid <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,w1,w2,w,AA1,AA2,Com,tol=1e-4,maxiter=500,report.prog=FALSE)
{

  n.lambda <- length(lambda.seq)
  n.eta <- length(eta.seq)

  b1.arr <- array(0,dim=c(ncol(X1),n.lambda,n.eta))
  b2.arr <- array(0,dim=c(ncol(X2),n.lambda,n.eta))
  
  P1.arr <- array(0,dim=c(nrow(X1),n.lambda,n.eta))
  P2.arr <- array(0,dim=c(nrow(X2),n.lambda,n.eta))
  
  iterations <- matrix(0,n.lambda*n.eta,3)
  colnames(iterations) <- c("lambda","eta","iter")
  step <- 0
  init <- list( beta1 = rep(0,ncol(X1)),
                beta2 = rep(0,ncol(X2)))
  
  
  for(l in 1:n.lambda){
    for(k in 1:n.eta){

      grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
                                                         rX1 = X1,
                                                         groups1 = groups1,
                                                         rY2 = Y2,
                                                         rX2 = X2,
                                                         groups2 = groups2,
                                                         rho1 = rho1,
                                                         rho2 = rho2,
                                                         lambda = lambda.seq[l],
                                                         eta = eta.seq[k],
                                                         w1 = w1,
                                                         w2 = w2,
                                                         w = w,
                                                         rAA1 = AA1,
                                                         rAA2 = AA2,
                                                         rCom = Com,
                                                         tol = tol,
                                                         maxiter = maxiter,
                                                         beta1_init = init$beta1,
                                                         beta2_init = init$beta2)

      b1 <- grouplasso2pop_logreg.out$beta1.hat
      b2 <- grouplasso2pop_logreg.out$beta2.hat

      init <- list( beta1 = b1,
                    beta2 = b2)

      b1.arr[,l,k] <- b1
      b2.arr[,l,k] <- b2

      P1.arr[,l,k] <- logit(X1 %*% b1)
      P2.arr[,l,k] <- logit(X2 %*% b2)
      
      step <- step + 1
      iterations[step,] <- c(lambda.seq[l],eta.seq[k],grouplasso2pop_logreg.out$iter)

      if(report.prog == TRUE){

        print(c(l,k,grouplasso2pop_logreg.out$iter))

      }

    }

  }

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  P1.arr = P1.arr,
                  P2.arr = P2.arr,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = iterations)

  return(output)

}


#' Choose tuning parameters by crossvalidation for grouplasso2pop logreg when given a fixed grid of lambda and eta values
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda.seq sequence of lambda values
#' @param eta.seq sequence of eta values
#' @param n.folds the number of crossvalidation folds
#' @param b1.init.arr array of initial values for beta1
#' @param b2.init.arr array of initial values for beta2
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @param tol the convergence tolerance
#' @param maxiter the maximum number of iterations allowed for each fit
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#' @examples
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 400,n2 = 600)
#' 
#' grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                              X1 = grouplasso2pop_logreg_data$X1,
#'                                                              groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                              Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                              X2 = grouplasso2pop_logreg_data$X2,
#'                                                              groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                              rho1 = 2,
#'                                                              rho2 = 1,
#'                                                              n.lambda = 5,
#'                                                              n.eta = 5,
#'                                                              lambda.min.ratio = 0.01,
#'                                                              w1 = grouplasso2pop_logreg_data$w1,
#'                                                              w2 = grouplasso2pop_logreg_data$w2,
#'                                                              w = grouplasso2pop_logreg_data$w,
#'                                                              AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                              AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                              Com = grouplasso2pop_logreg_data$Com,
#'                                                              tol = 1e-3,
#'                                                              maxiter = 500,
#'                                                              report.prog = TRUE)
#' 
#' lambda.seq <- grouplasso2pop_logreg_grid.out$lambda.seq
#' eta.seq <- grouplasso2pop_logreg_grid.out$eta.seq
#' b1.arr <- grouplasso2pop_logreg_grid.out$b1.arr
#' b2.arr <- grouplasso2pop_logreg_grid.out$b2.arr
#' 
#' grouplasso2pop_logreg_cv_fixedgrid.out <- grouplasso2pop_logreg_cv_fixedgrid(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                                              X1 = grouplasso2pop_logreg_data$X1,
#'                                                                              groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                                              Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                                              X2 = grouplasso2pop_logreg_data$X2,
#'                                                                              groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                                              rho1 = 2,
#'                                                                              rho2 = 1,
#'                                                                              lambda.seq = lambda.seq,
#'                                                                              eta.seq = eta.seq,
#'                                                                              n.folds = 5,
#'                                                                              b1.init.arr = b1.arr,
#'                                                                              b2.init.arr = b2.arr,
#'                                                                              w1 = grouplasso2pop_logreg_data$w1,
#'                                                                              w2 = grouplasso2pop_logreg_data$w2,
#'                                                                              w = grouplasso2pop_logreg_data$w,
#'                                                                              AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                                              AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                                              Com = grouplasso2pop_logreg_data$Com,
#'                                                                              tol = 1e-3,
#'                                                                              maxiter = 500)
#' @export
grouplasso2pop_logreg_cv_fixedgrid <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,n.folds,b1.init.arr,b2.init.arr,w1,w2,w,AA1,AA2,Com,tol=1e-3,maxiter=500)
{

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

        grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1[-fold1],
                                                           rX1 = X1[-fold1,],
                                                           groups1 = groups1,
                                                           rY2 = Y2[-fold2],
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
                                                           maxiter = maxiter,
                                                           beta1_init = b1.init.arr[,l,k],
                                                           beta2_init = b2.init.arr[,l,k])

        b1.fold <- grouplasso2pop_logreg.out$beta1.hat
        b2.fold <- grouplasso2pop_logreg.out$beta2.hat

        b1.folds.arr[,l,k,fold] <- b1.fold
        b2.folds.arr[,l,k,fold] <- b2.fold

        iterations[step,2+fold] <- grouplasso2pop_logreg.out$iter

        P1.fold <- logit(X1[fold1,] %*% b1.fold)
        P2.fold <- logit(X2[fold2,] %*% b2.fold)

        minus2ll1.fold <- - 2 * rho1 * mean( Y1[fold1] * log(P1.fold) + (1-Y1[fold1]) * log( 1 - P1.fold) )
        minus2ll2.fold <- - 2 * rho2 * mean( Y2[fold2] * log(P2.fold) + (1-Y2[fold2]) * log( 1 - P2.fold) )

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

#' Choose tuning parameters by crossvalidation for grouplasso2pop logreg.
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values desired
#' @param n.eta the number of eta values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#'
#' @examples
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 400,
#'                                                              n2 = 600)
#' 
#' grouplasso2pop_logreg_cv.out <- grouplasso2pop_logreg_cv(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                          X1 = grouplasso2pop_logreg_data$X1,
#'                                                          groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                          Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                          X2 = grouplasso2pop_logreg_data$X2,
#'                                                          groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                          rho1 = 2,
#'                                                          rho2 = 1,
#'                                                          n.lambda = 5,
#'                                                          n.eta = 5,
#'                                                          lambda.min.ratio = 0.01,
#'                                                          n.folds = 5,
#'                                                          w1 = grouplasso2pop_logreg_data$w1,
#'                                                          w2 = grouplasso2pop_logreg_data$w2,
#'                                                          w = grouplasso2pop_logreg_data$w,
#'                                                          AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                          AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                          Com = grouplasso2pop_logreg_data$Com,
#'                                                          tol = 1e-3,
#'                                                          maxiter = 500,
#'                                                          report.prog = TRUE)
#' @export
grouplasso2pop_logreg_cv <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,lambda.max.ratio=1,n.folds,w1,w2,w,AA1,AA2,Com,tol=1e-4,maxiter=500,report.prog = TRUE){

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set,
  # which will be used as initial values for the crossvalidation training fits.
  grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = Y1,
                                                               X1 = X1,
                                                               groups1 = groups1,
                                                               Y2 = Y2,
                                                               X2 = X2,
                                                               groups2 = groups2,
                                                               rho1 = rho1,
                                                               rho2 = rho2,
                                                               n.lambda = n.lambda,
                                                               n.eta = n.eta,
                                                               lambda.min.ratio = lambda.min.ratio,
                                                               lambda.max.ratio = lambda.max.ratio,
                                                               w1 = w1,
                                                               w2 = w2,
                                                               w = w,
                                                               AA1 = AA1,
                                                               AA2 = AA2,
                                                               Com = Com,
                                                               tol = tol,
                                                               maxiter = maxiter,
                                                               report.prog = report.prog)

  lambda.seq <- grouplasso2pop_logreg_grid.out$lambda.seq
  eta.seq <- grouplasso2pop_logreg_grid.out$eta.seq
  b1.arr <- grouplasso2pop_logreg_grid.out$b1.arr
  b2.arr <- grouplasso2pop_logreg_grid.out$b2.arr

  # do the crossvalidation
  grouplasso2pop_logreg_cv_fixedgrid.out <- grouplasso2pop_logreg_cv_fixedgrid(Y1 = Y1,
                                                                               X1 = X1,
                                                                               groups1 = groups1,
                                                                               Y2 = Y2,
                                                                               X2 = X2,
                                                                               groups2 = groups2,
                                                                               rho1 = rho1,
                                                                               rho2 = rho2,
                                                                               lambda.seq = lambda.seq,
                                                                               eta.seq = eta.seq,
                                                                               n.folds = n.folds,
                                                                               b1.init.arr = b1.arr,
                                                                               b2.init.arr = b2.arr,
                                                                               w1 = w1,
                                                                               w2 = w2,
                                                                               w = w,
                                                                               AA1 = AA1,
                                                                               AA2 = AA2,
                                                                               Com = Com,
                                                                               tol = tol,
                                                                               maxiter = 500)

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  b1.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplasso2pop_logreg_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  iterations = grouplasso2pop_logreg_cv_fixedgrid.out$iterations)

  return(output)

}


#' Choose tuning parameters by crossvalidation for grouplasso2pop logreg with adaptive weights
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param n.lambda the number of lambda values desired
#' @param n.eta the number of eta values desired
#' @param lambda.min.ratio ratio of the smallest lambda value to the smallest value of lambda which admits no variables to the model
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#'
#' @examples
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 500, n2 = 604)
#'
#' grouplasso2pop_logreg_cv_adapt.out <- grouplasso2pop_logreg_cv_adapt(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                                      X1 = grouplasso2pop_logreg_data$X1,
#'                                                                      groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                                      Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                                      X2 = grouplasso2pop_logreg_data$X2,
#'                                                                      groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                                      rho1 = 2,
#'                                                                      rho2 = 1,
#'                                                                      n.lambda = 5,
#'                                                                      n.eta = 5,
#'                                                                      lambda.min.ratio = 0.01,
#'                                                                      n.folds = 5,
#'                                                                      w1 = grouplasso2pop_logreg_data$w1,
#'                                                                      w2 = grouplasso2pop_logreg_data$w2,
#'                                                                      w = grouplasso2pop_logreg_data$w,
#'                                                                      AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                                      AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                                      Com = grouplasso2pop_logreg_data$Com,
#'                                                                      tol = 1e-2,
#'                                                                      maxiter = 500,
#'                                                                      report.prog = TRUE)
#' @export
grouplasso2pop_logreg_cv_adapt <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,n.lambda,n.eta,lambda.min.ratio,lambda.max.ratio=1,n.folds,w1,w2,w,AA1,AA2,Com,tol=1e-3,maxiter=500,report.prog = TRUE){

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
      
      norms1[j] <- sqrt(sum((t(X1[,ind1]) %*% (Y1 - mean(Y1)))^2)) / ( w1[j] * n1 / rho1)
      
    }
    
    if(j <= q2){
      
      norms2[j] <- sqrt(sum((t(X2[,ind2]) %*% (Y2 - mean(Y2)))^2)) / (w2[j] * n2 / rho1)
    }

  }

  # smallest value of lambda which sets all the non-intercept entries of beta1 and beta2 equal to zero.
  lambda.max <- 2 * max(norms1,norms2)
  lambda.initial.fit <- lambda.min.ratio * lambda.max

  # fit a grouplasso2pop with eta = 0 and lambda as lambda.min.ratio*lambda.max.
  
  grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
                                                     rX1 = X1,
                                                     groups1 = groups1,
                                                     rY2 = Y2,
                                                     rX2 = X2,
                                                     groups2 = groups2,
                                                     rho1 = rho1,
                                                     rho2 = rho2,
                                                     lambda = lambda.initial.fit,
                                                     eta = 0,
                                                     w1 = w1,
                                                     w2 = w2,
                                                     w = w,
                                                     rAA1 = AA1,
                                                     rAA2 = AA2,
                                                     rCom = Com,
                                                     tol = tol,
                                                     maxiter = maxiter)

  # now make new values of w1, w2, and w based on these.

  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    w1[j] <- min(w1[j]/sqrt( sum( grouplasso2pop_logreg.out$beta1.hat[ind]^2 )),1e10) #replace Inf with 1e10
  }
  for(j in 1:q2)
  {
    ind <- which(groups2 == j)
    w2[j] <- min(w2[j]/sqrt( sum( grouplasso2pop_logreg.out$beta2.hat[ind]^2 )),1e10)
  }
  for( j in Com)
  {
    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)
    w[j] <- min(w[j]/sum( (AA1[[j]] %*% grouplasso2pop_logreg.out$beta1.hat[ind1] - AA2[[j]] %*% grouplasso2pop_logreg.out$beta2.hat[ind2] )^2),1e10)
  }

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = Y1,
                                                               X1 = X1,
                                                               groups1 = groups1,
                                                               Y2 = Y2,
                                                               X2 = X2,
                                                               groups2 = groups2,
                                                               rho1 = rho1,
                                                               rho2 = rho2,
                                                               n.lambda = n.lambda,
                                                               n.eta = n.eta,
                                                               lambda.min.ratio = lambda.min.ratio,
                                                               lambda.max.ratio = lambda.max.ratio,
                                                               w1 = w1,
                                                               w2 = w2,
                                                               w = w,
                                                               AA1 = AA1,
                                                               AA2 = AA2,
                                                               Com = Com,
                                                               tol = tol,
                                                               maxiter = maxiter,
                                                               report.prog = TRUE)

  lambda.seq <- grouplasso2pop_logreg_grid.out$lambda.seq
  eta.seq <- grouplasso2pop_logreg_grid.out$eta.seq
  b1.arr <- grouplasso2pop_logreg_grid.out$b1.arr
  b2.arr <- grouplasso2pop_logreg_grid.out$b2.arr

  # do the crossvalidation
  grouplasso2pop_logreg_cv_fixedgrid.out <- grouplasso2pop_logreg_cv_fixedgrid(Y1 = Y1,
                                                                               X1 = X1,
                                                                               groups1 = groups1,
                                                                               Y2 = Y2,
                                                                               X2 = X2,
                                                                               groups2 = groups2,
                                                                               rho1 = rho1,
                                                                               rho2 = rho2,
                                                                               lambda.seq = lambda.seq,
                                                                               eta.seq = eta.seq,
                                                                               n.folds = n.folds,
                                                                               b1.init.arr = b1.arr,
                                                                               b2.init.arr = b2.arr,
                                                                               w1 = w1,
                                                                               w2 = w2,
                                                                               w = w,
                                                                               AA1 = AA1,
                                                                               AA2 = AA2,
                                                                               Com = Com,
                                                                               tol = tol,
                                                                               maxiter = maxiter)

  output <- list( b1.arr = b1.arr,
                  b2.arr = b2.arr,
                  P1.arr = grouplasso2pop_logreg_grid.out$P1.arr,
                  P2.arr = grouplasso2pop_logreg_grid.out$P2.arr,
                  b1.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplasso2pop_logreg_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w1 = w1,
                  w2 = w2,
                  w = w,
                  iterations = grouplasso2pop_logreg_cv_fixedgrid.out$iterations)

  return(output)

}

#' Choose tuning parameters by crossvalidation for grouplasso2pop logreg with adaptive weights
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 matrix containing the design matrices for data set 1
#' @param groups1 a vector indicating to which group each covariate of data set 1 belongs
#' @param Y2 the binary response vector of data set 2
#' @param X2 matrix containing the design matrices for data set 2
#' @param groups2 a vector indicating to which group each covariate of data set 2 belongs
#' @param rho1 weight placed on the first data set
#' @param rho2 weight placed on the second data set
#' @param lambda.seq the lambda sequence
#' @param eta.seq sequence of eta values
#' @param n.folds the number of crossvalidation folds
#' @param w1 group-specific weights for different penalization across groups in data set 1
#' @param w2 group-specific weights for different penalization across groups in data set 2
#' @param w group-specific weights for different penalization toward similarity for different groups
#' @param AA1 a list of the matrices A1j
#' @param AA1 a list of the matrices A2j
#' @param Com the indices of the covariate groups which are common in the two datasets
#' @return a list containing the fits over a grid of lambda and eta values as well as the vector of lambda values and the vector of eta values
#'
#' @examples
#' grouplasso2pop_logreg_data <- get_grouplasso2pop_logreg_data(n1 = 400, n2 = 600)
#' 
#' grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                              X1 = grouplasso2pop_logreg_data$X1,
#'                                                              groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                              Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                              X2 = grouplasso2pop_logreg_data$X2,
#'                                                              groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                              n.lambda = 5,
#'                                                              n.eta = 5,
#'                                                              lambda.min.ratio = 0.01,
#'                                                              w1 = grouplasso2pop_logreg_data$w1,
#'                                                              w2 = grouplasso2pop_logreg_data$w2,
#'                                                              w = grouplasso2pop_logreg_data$w,
#'                                                              AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                              AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                              Com = grouplasso2pop_logreg_data$Com,
#'                                                              tol = 1e-3,
#'                                                              maxiter = 500,
#'                                                              report.prog = TRUE)
#' 
#' grouplasso2pop_logreg_cv_adapt_fixedgrid.out <- grouplasso2pop_logreg_cv_adapt_fixedgrid(Y1 = grouplasso2pop_logreg_data$Y1,
#'                                                                                          X1 = grouplasso2pop_logreg_data$X1,
#'                                                                                          groups1 = grouplasso2pop_logreg_data$groups1,
#'                                                                                          Y2 = grouplasso2pop_logreg_data$Y2,
#'                                                                                          X2 = grouplasso2pop_logreg_data$X2,
#'                                                                                          groups2 = grouplasso2pop_logreg_data$groups2,
#'                                                                                          lambda.seq = grouplasso2pop_logreg_grid.out$lambda.seq,
#'                                                                                          eta.seq = grouplasso2pop_logreg_grid.out$eta.seq,
#'                                                                                          n.folds = 6,
#'                                                                                          lambda.initial.fit = grouplasso2pop_logreg_grid.out$lambda.seq[1],
#'                                                                                          w1 = grouplasso2pop_logreg_data$w1,
#'                                                                                          w2 = grouplasso2pop_logreg_data$w2,
#'                                                                                          w = grouplasso2pop_logreg_data$w,
#'                                                                                          AA1 = grouplasso2pop_logreg_data$AA1,
#'                                                                                          AA2 = grouplasso2pop_logreg_data$AA2,
#'                                                                                          Com = grouplasso2pop_logreg_data$Com)
#' @export
grouplasso2pop_logreg_cv_adapt_fixedgrid <- function(Y1,X1,groups1,Y2,X2,groups2,rho1,rho2,lambda.seq,eta.seq,lambda.initial.fit,n.folds,w1,w2,w,AA1,AA2,Com,tol=1e-3,maxiter=500,report.prog = TRUE)
{

  # get initial fit at small value of lambda and eta = 0

  grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
                                                     rX1 = X1,
                                                     groups1 = groups1,
                                                     rY2 = Y2,
                                                     rX2 = X2,
                                                     groups2 = groups2,
                                                     rho1 = rho1,
                                                     rho2 = rho2,
                                                     lambda = lambda.initial.fit,
                                                     eta = 0,
                                                     w1 = w1,
                                                     w2 = w2,
                                                     w = w,
                                                     rAA1 = AA1,
                                                     rAA2 = AA2,
                                                     rCom = Com,
                                                     tol = tol,
                                                     maxiter = maxiter)

  # now make new values of w1, w2, and w based on these.

  q1 <- length(unique(groups1))
  q2 <- length(unique(groups2))

  for(j in 1:q1)
  {
    ind <- which(groups1 == j)
    w1[j] <- min(w1[j]/sqrt( sum( grouplasso2pop_logreg.out$beta1.hat[ind]^2 )),1e10) #replace Inf with 1e10
  }
  for(j in 1:q2)
  {
    ind <- which(groups2 == j)
    w2[j] <- min(w2[j]/sqrt( sum( grouplasso2pop_logreg.out$beta2.hat[ind]^2 )),1e10)
  }
  for( j in Com)
  {
    ind1 <- which(groups1 == j)
    ind2 <- which(groups2 == j)
    w[j] <- min(w[j]/sum( (AA1[[j]] %*% grouplasso2pop_logreg.out$beta1.hat[ind1] - AA2[[j]] %*% grouplasso2pop_logreg.out$beta2.hat[ind2] )^2),1e10)
  }

  # obtain lambda.seq and eta.seq from the grid function, as well as the fits on the entire data set, which will be used as initial values for the crossvalidation training fits.
  grouplasso2pop_logreg_fixedgrid.out <- grouplasso2pop_logreg_fixedgrid(Y1 = Y1,
                                                                         X1 = X1,
                                                                         groups1 = groups1,
                                                                         Y2 = Y2,
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
                                                                         tol = tol,
                                                                         maxiter = maxiter,
                                                                         report.prog = report.prog)
                                                                    

  # do the crossvalidation
  grouplasso2pop_logreg_cv_fixedgrid.out <- grouplasso2pop_logreg_cv_fixedgrid(Y1 = Y1,
                                                                               X1 = X1,
                                                                               groups1 = groups1,
                                                                               Y2 = Y2,
                                                                               X2 = X2,
                                                                               groups2 = groups2,
                                                                               rho1 = rho1,
                                                                               rho2 = rho2,
                                                                               lambda.seq = lambda.seq,
                                                                               eta.seq = eta.seq,
                                                                               n.folds = n.folds,
                                                                               b1.init.arr = grouplasso2pop_logreg_fixedgrid.out$b1.arr,
                                                                               b2.init.arr = grouplasso2pop_logreg_fixedgrid.out$b2.arr,
                                                                               w1 = w1,
                                                                               w2 = w2,
                                                                               w = w,
                                                                               AA1 = AA1,
                                                                               AA2 = AA2,
                                                                               Com = Com,
                                                                               tol = tol,
                                                                               maxiter = maxiter)
  # collect output
  output <- list( b1.arr = grouplasso2pop_logreg_fixedgrid.out$b1.arr,
                  b2.arr = grouplasso2pop_logreg_fixedgrid.out$b2.arr,
                  P1.arr = grouplasso2pop_logreg_fixedgrid.out$P1.arr,
                  P2.arr = grouplasso2pop_logreg_fixedgrid.out$P2.arr,
                  b1.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b1.folds.arr,
                  b2.folds.arr = grouplasso2pop_logreg_cv_fixedgrid.out$b2.folds.arr,
                  minus2ll.arr = grouplasso2pop_logreg_cv_fixedgrid.out$minus2ll.arr,
                  which.lambda.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  lambda.seq = lambda.seq,
                  eta.seq = eta.seq,
                  lambda.initial.fit = lambda.initial.fit,
                  w1 = w1,
                  w2 = w2,
                  w = w,
                  iterations = grouplasso2pop_logreg_cv_fixedgrid.out$iterations)

  return(output)

}

#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity
#'
#' @param Y1 the binary response vector of data set 1
#' @param XX1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the binary response vector of data set 2
#' @param XX2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
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
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#'                                                          n2 = 604)
#'
#' semipadd2pop_logreg.out <- semipadd2pop_logreg(Y1 = semipadd2pop_logreg_data$Y1,
#'                                                X1 = semipadd2pop_logreg_data$X1,
#'                                                nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                Y2 = semipadd2pop_logreg_data$Y2,
#'                                                X2 = semipadd2pop_logreg_data$X2,
#'                                                nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                rho1 = 2,
#'                                                rho2 = 1,
#'                                                w1 = 1,
#'                                                w2 = 1,
#'                                                w = 1,
#'                                                nCom = semipadd2pop_logreg_data$nCom,
#'                                                d1 = semipadd2pop_logreg_data$nonparm1*25,
#'                                                d2 = semipadd2pop_logreg_data$nonparm2*15,
#'                                                xi = .5,
#'                                                lambda.beta = .01,
#'                                                lambda.f = .01,
#'                                                eta.beta = .01,
#'                                                eta.f = .01,
#'                                                tol = 1e-3,
#'                                                maxiter = 500,
#'                                                plot_obj = FALSE)
#'
#' plot_semipaddgt2pop(semipadd2pop_logreg.out,
#'                     true.functions=list(f1 = semipadd2pop_logreg_data$f1,
#'                                         f2 = semipadd2pop_logreg_data$f2,
#'                                         X1 = semipadd2pop_logreg_data$X1,
#'                                         X2 = semipadd2pop_logreg_data$X2)
#' )
#' @export
semipadd2pop_logreg <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,lambda.beta,lambda.f,eta.beta,eta.f,tol=1e-4,maxiter=500,plot_obj=FALSE)
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
  grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = Y1,
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

  # construct fitted functions from grouplasso2pop output
  semipadd2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
                                                        nonparm1 = nonparm1,
                                                        groups1 = grouplasso2pop_inputs$groups1,
                                                        knots.list1 = grouplasso2pop_inputs$knots.list1,
                                                        emp.cent1 = grouplasso2pop_inputs$emp.cent1,
                                                        QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
                                                        b1 = grouplasso2pop_logreg.out$beta1.hat,
                                                        X2 = X2,
                                                        nonparm2 = nonparm2,
                                                        groups2 = grouplasso2pop_inputs$groups2,
                                                        knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                        emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                        QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                        b2 = grouplasso2pop_logreg.out$beta2.hat)

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

  class(output) <- "semipaddgt2pop"

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
semipadd2pop_logreg_grid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
  grouplasso2pop_logreg_grid.out <- grouplasso2pop_logreg_grid(Y1 = Y1,
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
                                                               w1 = grouplasso2pop_inputs$w1,
                                                               w2 = grouplasso2pop_inputs$w2,
                                                               w = grouplasso2pop_inputs$w,
                                                               AA1 = grouplasso2pop_inputs$AA1.tilde,
                                                               AA2 = grouplasso2pop_inputs$AA2.tilde,
                                                               Com = grouplasso2pop_inputs$Com,
                                                               tol = tol,
                                                               maxiter = maxiter,
                                                               report.prog = report.prog)

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
                                                              b1 = grouplasso2pop_logreg_grid.out$b1[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_logreg_grid.out$b2[,l,k])

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
                  P1.hat = grouplasso2pop_logreg_grid.out$P1.arr,
                  P2.hat = grouplasso2pop_logreg_grid.out$P2.arr,
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
                  lambda.seq = grouplasso2pop_logreg_grid.out$lambda.seq,
                  eta.seq = grouplasso2pop_logreg_grid.out$eta.seq,
                  iterations = grouplasso2pop_logreg_grid.out$iterations)

  class(output) <- "semipaddgt2pop_grid"

  return(output)

}



#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the binary response vector of data set 2
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
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#'                                                          n2 = 604)
#'
#' semipadd2pop_logreg_cv.out <- semipadd2pop_logreg_cv(Y1 = semipadd2pop_logreg_data$Y1,
#'                                                      X1 = semipadd2pop_logreg_data$X1,
#'                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                      Y2 = semipadd2pop_logreg_data$Y2,
#'                                                      X2 = semipadd2pop_logreg_data$X2,
#'                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                      rho1 = 2,
#'                                                      rho2 = 1,
#'                                                      w1 = 1,
#'                                                      w2 = 1,
#'                                                      w = 1,
#'                                                      nCom = semipadd2pop_logreg_data$nCom,
#'                                                      d1 = semipadd2pop_logreg_data$nonparm1*25,
#'                                                      d2 = semipadd2pop_logreg_data$nonparm2*15,
#'                                                      xi = .5,
#'                                                      n.lambda = 5,
#'                                                      n.eta = 5,
#'                                                      lambda.min.ratio = .01,
#'                                                      n.folds = 5,
#'                                                      lambda.beta = 1,
#'                                                      lambda.f = 1,
#'                                                      eta.beta = 1,
#'                                                      eta.f = 1,
#'                                                      tol = 1e-3,
#'                                                      maxiter = 1000,
#'                                                      report.prog = FALSE)
#'
#' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv.out,
#'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#'                                              f2 = semipadd2pop_logreg_data$f2,
#'                                              X1 = semipadd2pop_logreg_data$X1,
#'                                              X2 = semipadd2pop_logreg_data$X2)
#' )
#' @export
semipadd2pop_logreg_cv <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1,w2,w,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio=.01,lambda.max.ratio=1,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
  grouplasso2pop_logreg_cv.out <- grouplasso2pop_logreg_cv(Y1 = Y1,
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
                                                              b1 = grouplasso2pop_logreg_cv.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_logreg_cv.out$b2.arr[,l,k])

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
                                                                b1 = grouplasso2pop_logreg_cv.out$b1.folds.arr[,l,k,fold],
                                                                X2 = X2,
                                                                nonparm2 = nonparm2,
                                                                groups2 = grouplasso2pop_inputs$groups2,
                                                                knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                                emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                                QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                                b2 = grouplasso2pop_logreg_cv.out$b2.folds.arr[,l,k,fold])

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
                  lambda.seq = grouplasso2pop_logreg_cv.out$lambda.seq,
                  eta.seq = grouplasso2pop_logreg_cv.out$eta.seq,
                  which.lambda.cv = grouplasso2pop_logreg_cv.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv.out$which.lambda.cv.under.zero.eta,
                  iterations = grouplasso2pop_logreg_cv.out$iterations)

  class(output) <- "sempaddgt2pop_cv"

  return(output)

}

#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the binary response vector of data set 2
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
#' @param lambda.max.ratio ratio of the largest lambda value to the smallest value of lambda which admits no variables to the model
#' @param n.folds the number of crossvalidation folds
#' @param lambda.beta the level of sparsity penalization for the parametric effects (relative to nonparametric effects)
#' @param lambda.f the level of sparsity penalization for the nonparametric effects (relative to the parametric effects)
#' @param eta.beta the level of penalization towards model similarity for parametric effects indicated to be common (relative to nonparametric effects)
#' @param eta.f the level of penalization towards model similarity for nonparametric effects indicated to be common (relative to the parametric effects)
#' @param tol a convergence criterion
#' @param maxiter the maximum allowed number of iterations
#' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
#' @return Returns the estimator of the semiparametric additive model
#'
#' @examples
#' semipadd2pop_logreg_data <- get_semipadd2pop_logreg_data(n1 = 501,
#'                                                          n2 = 604)
#'
#' semipadd2pop_logreg_cv.out <- semipadd2pop_logreg_cv(Y1 = semipadd2pop_logreg_data$Y1,
#'                                                      X1 = semipadd2pop_logreg_data$X1,
#'                                                      nonparm1 = semipadd2pop_logreg_data$nonparm1,
#'                                                      Y2 = semipadd2pop_logreg_data$Y2,
#'                                                      X2 = semipadd2pop_logreg_data$X2,
#'                                                      nonparm2 = semipadd2pop_logreg_data$nonparm2,
#'                                                      rho1 = 2,
#'                                                      rho2 = 1,
#'                                                      w1 = 1,
#'                                                      w2 = 1,
#'                                                      w = 1,
#'                                                      nCom = semipadd2pop_logreg_data$nCom,
#'                                                      d1 = semipadd2pop_logreg_data$nonparm1*25,
#'                                                      d2 = semipadd2pop_logreg_data$nonparm2*15,
#'                                                      xi = .5,
#'                                                      n.lambda = 5,
#'                                                      n.eta = 5,
#'                                                      lambda.min.ratio = .001,
#'                                                      n.folds = 5,
#'                                                      lambda.beta = 1,
#'                                                      lambda.f = 1,
#'                                                      eta.beta = 1,
#'                                                      eta.f = 1,
#'                                                      tol = 1e-3,
#'                                                      maxiter = 1000,
#'                                                      report.prog = FALSE)
#'
#' plot_semipaddgt2pop_cv(semipadd2pop_logreg_cv.out,
#'                        true.functions = list(f1 = semipadd2pop_logreg_data$f1,
#'                                              f2 = semipadd2pop_logreg_data$f2,
#'                                              X1 = semipadd2pop_logreg_data$X1,
#'                                              X2 = semipadd2pop_logreg_data$X2)
#' )
#' @export
semipadd2pop_logreg_cv_adapt <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,w1=1,w2=1,w=1,nCom,d1,d2,xi,n.lambda = 5,n.eta = 5,lambda.min.ratio = .01,lambda.max.ratio=1,n.folds = 5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
  grouplasso2pop_logreg_cv_adapt.out <- grouplasso2pop_logreg_cv_adapt(Y1 = Y1,
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
                                                              b1 = grouplasso2pop_logreg_cv_adapt.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2 = nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_logreg_cv_adapt.out$b2.arr[,l,k])

      f1.hat[[l]][[k]] <- semipaddgt2pop_fitted$f1.hat
      f2.hat[[l]][[k]] <- semipaddgt2pop_fitted$f2.hat

      f1.hat.design[,,l,k] <- semipaddgt2pop_fitted$f1.hat.design
      f2.hat.design[,,l,k] <- semipaddgt2pop_fitted$f2.hat.design

      beta1.hat[,l,k] <- semipaddgt2pop_fitted$beta1.hat
      beta2.hat[,l,k] <- semipaddgt2pop_fitted$beta2.hat

      # for( fold in 1:n.folds)
      # {
      # 
      #   semipaddgt2pop_fitted <- grouplasso2pop_to_semipadd2pop(X1 = X1,
      #                                                           nonparm1 = nonparm1,
      #                                                           groups1 = grouplasso2pop_inputs$groups1,
      #                                                           knots.list1 = grouplasso2pop_inputs$knots.list1,
      #                                                           emp.cent1 = grouplasso2pop_inputs$emp.cent1,
      #                                                           QQ1.inv = grouplasso2pop_inputs$QQ1.inv,
      #                                                           b1 = grouplasso2pop_logreg_cv_adapt.out$b1.folds.arr[,l,k,fold],
      #                                                           X2 = X2,
      #                                                           nonparm2 = nonparm2,
      #                                                           groups2 = grouplasso2pop_inputs$groups2,
      #                                                           knots.list2 = grouplasso2pop_inputs$knots.list2,
      #                                                           emp.cent2 = grouplasso2pop_inputs$emp.cent2,
      #                                                           QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
      #                                                           b2 = grouplasso2pop_logreg_cv_adapt.out$b2.folds.arr[,l,k,fold])
      # 
      #   f1.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f1.hat
      #   f2.hat.folds[[l]][[k]][[fold]] <- semipaddgt2pop_fitted$f2.hat
      # 
      # }

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
                  P1.hat = grouplasso2pop_logreg_cv_adapt.out$P1.arr,
                  P2.hat = grouplasso2pop_logreg_cv_adapt.out$P2.arr,
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
                  lambda.seq = grouplasso2pop_logreg_cv_adapt.out$lambda.seq,
                  eta.seq = grouplasso2pop_logreg_cv_adapt.out$eta.seq,
                  lambda.initial.fit = grouplasso2pop_logreg_cv_adapt.out$lambda.initial.fit,
                  which.lambda.cv = grouplasso2pop_logreg_cv_adapt.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv_adapt.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_adapt.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplasso2pop_logreg_cv_adapt.out$w1,
                  w2 = grouplasso2pop_logreg_cv_adapt.out$w2,
                  w = grouplasso2pop_logreg_cv_adapt.out$w,
                  iterations = grouplasso2pop_logreg_cv_adapt.out$iterations)

  class(output) <- "semipaddgt_cv"

  return(output)

}

#' Compute semiparametric binary-response regression model with 2 data sets while penalizing dissimilarity using CV to select tuning parameters after an adaptive step
#'
#' @param Y1 the binary response vector of data set 1
#' @param X1 the matrix with the observed covariate values for data set 1 (including a column of ones for the intercept)
#' @param nonparm1 a vector indicating for which covariates a nonparametric function is to be estimated for data set 1
#' @param Y2 the binary response vector of data set 2
#' @param X2 the matrix with the observed covariate values for data set 2 (including a column of ones for the intercept)
#' @param nonparm2 a vector indicating for which covariates a nonparametric function is to be estimated for data set 2
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
semipadd2pop_logreg_cv_adapt_fixedgrid <- function(Y1,X1,nonparm1,Y2,X2,nonparm2,rho1,rho2,nCom,d1,d2,xi,w1,w2,w,lambda.seq,eta.seq,lambda.initial.fit,n.folds=5,lambda.beta=1,lambda.f=1,eta.beta=1,eta.f=1,tol=1e-3,maxiter = 1000,report.prog = FALSE)
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
  grouplasso2pop_logreg_cv_adapt_fixedgrid.out <- grouplasso2pop_logreg_cv_adapt_fixedgrid(Y1 = Y1,
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
  
  which.lambda.cv <- grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv
  which.eta.cv <- grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.eta.cv
  
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
                                                              b1 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$b1.arr[,l,k],
                                                              X2 = X2,
                                                              nonparm2,
                                                              groups2 = grouplasso2pop_inputs$groups2,
                                                              knots.list2 = grouplasso2pop_inputs$knots.list2,
                                                              emp.cent2 = grouplasso2pop_inputs$emp.cent2,
                                                              QQ2.inv = grouplasso2pop_inputs$QQ2.inv,
                                                              b2 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$b2.arr[,l,k])

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
                  P1.hat = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$P1.arr,
                  P2.hat = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$P2.arr,
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
                  which.lambda.cv = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv,
                  which.eta.cv = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.eta.cv,
                  which.lambda.cv.under.zero.eta = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$which.lambda.cv.under.zero.eta,
                  w1 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w1,
                  w2 = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w2,
                  w = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$w,
                  iterations = grouplasso2pop_logreg_cv_adapt_fixedgrid.out$iterations)

  class(output) <- "semipadd2pop_gt_cv"

  return(output)

}
