#' @title Get design matrix for finding regression splines coefficients under wiggliness penalization
#' 
#' @param X a vector of values at which penalized cubic B-splines are to be evaluated
#' @param d the number of B-spline functions in the basis.  If negative, knots are evenly spaced.
#' @param xi a parameter for penalizing towards smoothness
#' @param W a vector of observation weights
#' @return a list with the knots, the empirical centerings of the columns, and some matrices
#' @export
spsm_cubespline_design <- function(X,d,xi,W = NULL)
{
  
  n <- length(X)
  
  if(d < 0){
    
    int.knots <- seq(min(X), max(X), length = - d - 2 + 1 )
    d <- - d
    
  } else {
    
    int.knots <- int.knots.at.empirical.quantiles(X,d)
    # int.knots <- quantile(X,seq(0,1,length = d - 2 + 1)) # add one, so that one can be removed after centering to restore full-rank.
    
  }
  
  boundary.knots <- range(int.knots)
  all.knots <- sort(c(rep(boundary.knots,3),int.knots))
  
  if( length(W) == 0){
    
    B <- spline.des(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)$design[,-1] # remove one so we can center and keep full-rank
    emp.cent <- apply(B,2,mean)
    
  } else { # W is a vector of weights for each observation (has 0s and 1s in our context)
    
    B <- diag(W) %*% splineDesign(all.knots,X,ord=4,derivs=0,outer.ok=TRUE)[,-1] # remove one so we can center and keep full-rank
    emp.cent <- apply(B[which(W!=0), ],2,mean)
  }
  
  B.cent <- B - matrix(emp.cent,n,d,byrow=TRUE)
  
  # construct matrix in which l2 norm of function is a quadratic form
  M <- t(B.cent) %*% B.cent / n
  
  # construct matrix in which 2nd derivative penalty is a quadratic form
  # R <- matrix(NA,d+1,d+1)
  # dsq_bspline.mat <- splineDesign(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)
  # for(k in 1:(d+1))
  #   for(l in 1:(d+1))
  #   {
  # 
  #     pcwiselin <- dsq_bspline.mat[,k] * dsq_bspline.mat[,l] # Get sum of trapezoidal areas.
  #     h <- diff(int.knots)
  #     R[k,l] <- sum(.5*(pcwiselin[-1] + pcwiselin[-length(int.knots)])*h)  # sum of trapezoidal areas.
  # 
  #   }
  
  # # construct matrix in which 2nd derivative penalty is a quadratic form
  R <- matrix(NA,d+1,d+1)
  dsq_bspline.mat <- splineDesign(int.knots,knots = all.knots,outer.ok=TRUE,derivs=2)
  K <- length(int.knots)
  for(l in 1:(d+1))
    for(j in 1:(d+1)){

      tk <- dsq_bspline.mat[-K,l]
      gk <- dsq_bspline.mat[-K,j]

      tkp1 <- dsq_bspline.mat[-1,l]
      gkp1 <- dsq_bspline.mat[-1,j]

      R[l,j] <- sum( diff(int.knots) * ((1/2) * (gk*tkp1 + gkp1*tk) + (1/3)* (tkp1 - tk)*( gkp1 - gk) ))

    }
  
  
  Q <- chol(M + xi^2 * R[-1,-1]) # remove the one corresponding to the first coefficient (since we have removed one column)
  
  Q.inv <- solve(Q)
  D.tilde <- B.cent %*% Q.inv
  
  output <- list(knots = all.knots,
                 emp.cent = emp.cent,
                 Q.inv = Q.inv,
                 D.tilde = D.tilde,
                 d = d)
  
}

#' Put knots at empirical quantiles, accommodating data with ties.
#' 
#' @param X the data
#' @param d the number of basis functions desired (cubic B-splines)
#' @return the interior knots needed for constructing a centered, full-rank, cubic B-spline basis with \code{d} functions
#' @export
int.knots.at.empirical.quantiles <-  function(X,d){
  
  d_star <- d - 2 + 1
  d_try <- d_star
  n.int.knots <- 0
  
  while( n.int.knots < d_star)
  {
    
    
    int.knots <- unique(quantile(X,seq(0,1,length = d_try)))
    n.int.knots <- length(unique(int.knots))
    d_try <- d_try + 1
    
  }
  
  # if you have gotten too many knots, which will happen only rarely, remove some at random.
  if(n.int.knots > d_star){
    
    knots.rm <- sample(2:(n.int.knots - 1),n.int.knots - d_star,replace = FALSE)
    int.knots <- int.knots[-knots.rm]
    
  }
  
  return(int.knots)
  
}

#' The logit function
#'
#' @param z The value at which the logit function is to be evaluated.
#' @return the value of the logit function
#' @export
logit <- function(z){
  return(as.numeric(1/(1+exp(-z))))
}




#' Generate correlated Uniform(0,1) random variables
#'
#' @param n the desired number of sets of correlated Uniform random variables to generate
#' @param Rho the desired correlation matrix
#' @return a matrix containing the realizations in the rows
#'
#' @examples
#' n <- 1000
#' q <- 10
#' X <- corrUnif(n,Rho = (10/20)^abs(outer(1:q,1:q,"-")))
#'
#' Rho
#' cor(X)
#' @export
corrUnif <- function(n,Rho)
{
  q <- nrow(Rho)
  P <- 2*sin( Rho * pi / 6)
  X <- pnorm( matrix(rnorm(n* q),ncol= q) %*% chol(P))

  return(X)

}

#' Generate Bernoulli random vector with a certain correlation matrix
#'
#' @param n the desired number of Bernoulli random vectors to generate
#' @param probs the marginal probabilities of success for the Bernoulli random variables
#' @param Rho the desired correlation matrix (error will occur if the correlation matrix is inadmissable).
#' @return a matrix of which the rows are the Bernoulli random vectors
#'
#' This function implements the method from Chul Gyu Park, Taesung Park and Dong Wan Shin, A Simple Method for Generating Correlated Binary Variates
#' The American Statistician, Vol. 50, No. 4 (Nov., 1996), pp. 306-310, for generating correlated Bernoulli random variables.
#'
#' @examples
#' p <- 6 # number of Bernoulli random variables
#' Rho <- .2^abs(outer(1:p,1:p,"-"))  # matrix of desired correlations (must all be non-negative!)
#' probs <- 1:p/(2*p) # vector of marginal probabilities
#' n <- 500 # the number of independent realizations of the Bernoulli random vector
#'
#' Z <- corrBern(n,probs,Rho)
#'
#' Rho
#' cor(Z)
#' @export
corrBern <- function(n,probs,Rho)
{

  p <- ncol(Rho)
  A <- matrix(0,p,p)
  for(i in 1:p)
    for(j in i:p)
    {

      if(i == j)
      {

        A[i,i] <- -log(probs[i])

      } else if(i != j)
      {

        A[i,j] <- log( Rho[i,j] * sqrt( (1-probs[i])/probs[i]*(1-probs[j])/probs[j] ) + 1 )
        A[j,i] <- A[i,j]

      }

    }

  # check feasibility of correlation matrix under given probs:

  a <- numeric()
  for( i in 1:p)
  {
    a[i] <- 2*A[i,i] - sum(A[i,])
  }

  if(any(a < 0))
  {
    stop("The the given correlation matrix and marginal probabilities are incompatible")
  }


  X <- array(0,dim=c(p,p,n))
  for(i in 1:p)
    for(j in i:p)
    {

      if(i == j){

        X[i,i,] <- rpois(n,2*A[i,i] - sum(A[i,]))

      } else if(i != j)
      {

        X[i,j,] <- rpois(n,A[i,j])
        X[j,i,] <- X[i,j,]

      }

    }

  Y <- matrix(NA,n,p)
  for(k in 1:n)
    for(i in 1:p)
    {

      Y[k,i] <- sum(X[i,,k])

    }

  mean(Y[,1]==0)
  mean(Y[,2]==0)

  Z <- 1*(Y == 0)

  return(Z)

}


#' Generate a data set for group lasso
#'
#' @param n the sample size
#' @return a list containing the data
#' @export
get_grouplasso_data <- function(n,response){
  
  d <- c(1,1,1,1,1,3,4,10,20)
  q <- length(d)
  X <- matrix(rnorm(n*sum(d)),n,sum(d))
  groups <- numeric() ; for(j in 1:q){ groups <- c(groups,rep(j,d[j])) }
  beta <- c(0,2,0,0,0,1,1,1,0,0,0,0,rep(-1,10),rep(0,20))
  
  # set tuning parameters
  w <- rexp(q,2)
  
  if( response == "continuous"){
    
    Y <- X %*% beta + rnorm(n)
    
  } else if(response == "binary"){
    
    P <- logit(X %*% beta)
    Y <- rbinom(n,1,P)
   
  } else if(response == "gt"){
    
    # generate true response values
    P <- logit(X %*% beta)
    Y.true <- rbinom(n,1,P)
    
    # generate dorfman testing outcomes
    Se <- c(.98,.96)
    Sp <- c(.97,.99)
    assay.out <- dorfman.assay.gen(Y.true,Se,Sp,cj=4)
    Y <- list( A = assay.out$Z,
                I = assay.out$Y,
                Se = Se,
                Sp = Sp,
                E.approx = FALSE)
  
  } else {
    
    stop("invalid response type")
    
  }
  
  output <- list(Y = Y,
                 X = X,
                 groups = groups,
                 w = w,
                 beta = beta)
}



#' Generate two data sets with some common covariates
#'
#' @param n1 the sample size for the first data set
#' @param n2 the sample size for the second data set
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @return a list containing the data
#' @export
get_grouplasso2pop_data <- function(n1,n2,response){
  
  d1 <- c(1,1,1,4,5,rep(4,20))
  q1 <- length(d1)
  X1 <- matrix(rnorm(n1*sum(d1)),n1,sum(d1))
  groups1 <- numeric() ; for(j in 1:q1){ groups1 <- c(groups1,rep(j,d1[j])) }
  beta1 <- rnorm(ncol(X1))
  
  d2 <- c(1,1,4,3,4,1,rep(4,20))
  q2 <- length(d2)
  X2 <- matrix(rnorm(n2*sum(d2)),n2,sum(d2))
  groups2 <- numeric() ; for(j in 1:q2){ groups2 <- c(groups2,rep(j,d2[j])) }
  beta2 <- rnorm(ncol(X2))
  
  
  # set tuning parameters
  Com <- c(3,4)
  AA1 <- AA2 <- vector("list",min(q1,q2))
  for(j in Com){
    n.int <- rpois(1,4)
    AA1[[j]] <- matrix(rnorm(n.int*d1[j]),n.int,d1[j])
    AA2[[j]] <- matrix(rnorm(n.int*d2[j]),n.int,d2[j])
  }
  
  w1 <- rexp(q1,2)
  w2 <- rexp(q2,2)
  w <- rexp(min(q1,q2),2)
  
  
  if( response == "continuous"){
    

    Y1 <- X1 %*% beta1 + rnorm(n1)
    Y2 <- X2 %*% beta2 + rnorm(n2)
  
  } else if(response == "binary"){
    
    
    P1 <- logit(X1 %*% beta1)
    Y1 <- rbinom(n1,1,P1)
    
    P2 <- logit(X2 %*% beta2)
    Y2 <- rbinom(n2,1,P2)
    
    
  } else if(response == "gt"){
    
    
    # generate true response values
    P1 <- logit(X1 %*% beta1)
    Y1.true <- rbinom(n1,1,P1)
    
    # generate dorfman testing outcomes
    Se1 <- c(.98,.96)
    Sp1 <- c(.97,.99)
    assay1.out <- dorfman.assay.gen(Y1.true,Se1,Sp1,cj=4)
    Y1 <- list( A = assay1.out$Z,
                I = assay1.out$Y,
                Se = Se1,
                Sp = Sp1,
                E.approx = FALSE)
    
    # generate true response values
    P2 <- logit(X2 %*% beta2)
    Y2.true <- rbinom(n2,1,P2)
    
    # generate individual testing outcomes
    Se2 <- .96
    Sp2 <- .99
    assay2.out <- individual.assay.gen(Y2.true,Se2,Sp2,cj=1)
    Y2 <- list( A = assay2.out$Z,
                I = assay2.out$Y,
                Se = Se2,
                Sp = Sp2,
                E.approx = FALSE)
  
  } else {
    
    stop("invalid response type")
    
  }
  
  
  output <- list(Y1 = Y1,
                 X1 = X1,
                 groups1 = groups1,
                 Y2 = Y2,
                 X2 = X2,
                 groups2 = groups2,
                 w1 = w1,
                 w2 = w2,
                 w = w,
                 AA1 = AA1,
                 AA2 = AA2,
                 Com = Com,
                 response = response)
  
}




#'Generates individual testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se The testing sensitivity used for individual testing.
#' @param Sp The testing specificity used for individual testing.
#' @param cj This is an extraneous argument and is here only in order that all four assay.gen functions take the same arguments. Default is \code{cj=1} and an error will be returns if \code{cj!=1}.
#' @return a list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates individual level testing and stores the
#' testing responses in accordance to the data structure
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure
#' see McMahan et al. (2017).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' Y.true <- rbinom(100,1,p=.05)
#' assay.data <- individual.assay.gen(Y.true,Se=.96,Sp=.98)
#' @export
individual.assay.gen <-function(Y.true,Se,Sp,cj=1){
  if(cj!=1 || max(length(Se),length(Sp))>1)
  {
    print("cj must be 1 for individual testing and Se and Sp must each be of length 1")

  } else {

    N <- length(Y.true)

    Y<-matrix(-99,nrow=N, ncol=3) 		# initialize output matrix Y
    Z<-matrix(-99,nrow=N,ncol=cj+3) 	# initialize output matrix Z

    Z[,1] <- rbinom(N,1,(Se*Y.true+(1-Sp)*(1-Y.true)))
    Z[,2] <- 1
    Z[,3] <- 1
    Z[,4] <- 1:N

    Y[,2] <- 1
    Y[,3] <- 1:N
    return(list("Z"=Z, "Y"=Y))
  }

}

#'Generates master pool testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se The master pool testing sensitivity.
#' @param Sp The master pool testing specificity.
#' @param cj The size of the master pools (Note: The number of individuals \code{length(Y.true)} should be
#'      evenly divisible by \code{cj}, this is only for decoding purposes; i.e.,
#'      the regression methods do not require this condition).
#' @return a list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates Initial pool testing and stores the
#' testing responses in accordance to the data structure
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure
#' see McMahan et al. (2017).
#'
#' @examples
#' Y.true <- rbinom(100,1,p=.05)
#' assay.data <- masterpool.assay.gen(Y.true,Se=.96,Sp=.98,cj=4)
#' @export
masterpool.assay.gen <-function(Y.true,Se,Sp,cj){

  if(length(Se)>1 | length(Sp)>1)
  {

    print("error: Se or Sp has length greater than 1")
    return(NA)

  } else {

    N<-length(Y.true)
    Jmax<-N/cj
    J<-1

    Y<-matrix(-99,nrow=N, ncol=3)
    Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


    for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
    }

    J<-J-1
    Z<-Z[1:J,]

    return(list("Z"=Z, "Y"=Y))

  }
}

#'Generates Dorfman testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing sensitivity for the master pools and the second entry is the
#'      test sensitivity for individual testing.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for the master pools and the second entry is the
#'      test specificity for individual testing.
#' @param cj The size of the master pools (Note: The number of individuals \code{length(Y.true)} should be
#'      evenly divisible by \code{cj}. This is only for decoding purposes; i.e.,
#'      the regression methods do not require this condition).
#' @return A list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates Dorfman decoding and stores the
#' testing responses in accordance to the data structure
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure
#' see McMahan et al. (2017).
#'
#' @examples
#' Y.true <- rbinom(100,1,p=.05)
#' assay.data <- dorfman.assay.gen(Y.true,Se=c(.95,.92),Sp=c(.97,.98),cj=4)
#' @export
dorfman.assay.gen <-function(Y.true,Se,Sp,cj){
  N<-length(Y.true)
  Jmax<-N+N/cj
  J<-1

  Y<-matrix(-99,nrow=N, ncol=4)
  Z<-matrix(-99,nrow=Jmax,ncol=cj+3)
  D<-numeric(N)


  for(j in 1:(N/cj)){
    prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
    Z[J,1]<-rbinom(n=1,size=1,prob=prob)
    Z[J,2]<-cj
    Z[J,3]<-1
    Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
    Y[((j-1)*cj+1):(j*cj),2]<-1
    Y[((j-1)*cj+1):(j*cj),3]<-J
    J<-J+1
    if(Z[J-1,1]==1){
      for(k in ((j-1)*cj+1):(j*cj)){
        prob<-ifelse(Y.true[k]>0,Se[2],1-Sp[2])
        Z[J,1]<- rbinom(n=1,size=1,prob=prob)
        Z[J,2]<-1
        Z[J,3]<-2
        Z[J,4]<-k
        Y[k,2]<-2
        Y[k,4]<-J
        D[k] <- Z[J,1] # store individual assay results
        J<-J+1
      }
    } else { D[((j-1)*cj+1):(j*cj)]<-0 } # store individual assay results
  }

  J<-J-1
  Z<-Z[1:J,]

  return(list("Z"=Z, "Y"=Y, "D" = D))
}

#'Generates array testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing sensitivity for the row/column pools and the second entry is the
#'      test sensitivity for individual testing.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for the row/column pools and the second entry is the
#'      test specificity for individual testing.
#' @param cj Row and column pool sizes to be used (Note: The number of individuals
#'      should be evenly divisible by \code{cj*cj}. This is only for decoding
#'      purposes; i.e., the regression methods do not require this condition)
#' @return A list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates array decoding and stores the
#' testing responses in accordance to the data structure
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure
#' see McMahan et al. (2017).
#'
#' @examples
#' Y.true <- rbinom(100,1,p=.05)
#' assay.data <- array.assay.gen(Y.true,Se=c(.95,.92),Sp=c(.97,.98),cj=5)
#' @export
array.assay.gen <- function(Y.true, Se, Sp, cj){
  N<-length(Y.true) 					# get number of individuals
  Jmax<-2*N/cj + N					# specify maximum number of assays
  J<-1							# initialize assay index
  AT<-N/(cj^2)						# number of arrays (so N must be divisible by cj^2)

  Y<-matrix(-99,nrow=N, ncol=5) 			# initialize output matrix Y
  Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 			# initialize output matrix Z

  Y.A<-array(-99,c(cj,cj,AT))				# initialize array object to contain true individual statuses in the arrays
  ID.A<-array(-99,c(cj,cj,AT))				# initialize array object to contain individual indices in the array tests
  ind<-1
  for(i in 1:AT){
    for(j in 1:cj){
      for(m in 1:cj){
        Y.A[m,j,i]<-Y.true[ind]				# populate Y.A with the true statuses for the individuals in each array
        ID.A[m,j,i]<-ind					# populate ID.A with the IDs of the individuals in each array
        ind<-ind+1
      }}}

  D <- numeric(N)					#*** initialize vector to be populated with the individual diagnoses

  for(s in 1:AT){					# begin loop through arrays

    array.yk<-Y.A[,,s]					# pull the true statuses from array s
    array.id<-ID.A[,,s]					# pull the individual IDs from array s

    a<-rep(0,nrow(array.yk))				# initialize vector in which to store the row assays on array s
    b<-rep(0,ncol(array.yk))				# initialize vector in which to store the column assays on array s

    for(i in 1:cj){					# carry out row assays on array s
      prob<- ifelse(sum(array.yk[i,])>0, Se[1], 1-Sp[1]) # compute probability of a positive assay on row i of array s
      g<- rbinom(n=1,size=1,prob=prob)			# generate random assay on row i of array s
      a[i]<-g						# store randomly generated assay on row i of array s in position i of the vector a
      Z[J,1]<-g 						# store this also in row J of Z, column 1
      Z[J,2]<-cj 						# store the pool size of this assay in row J of Z, column 2
      Z[J,3]<-1						# store a 1 in row J of Z, column 3, to be used later to reference Se[1], Sp[1]
      Z[J,4:(cj+3)]<-array.id[i,]			# store the indices of the individuals in row i of array s in row J of Z, in the remaining cj columns
      Y[array.id[i,],2]<-2				# store a 2 in the rows of Y corresponding to the individuals in row i of array s; they are tested twice (once in a row, once in a column)
      Y[array.id[i,],3]<-J				# store J, the index of the assay, in the rows of Y corresponding to the individuals in row i of array s
      J<-J+1						# increment assay index
    }
    for(j in 1:cj){					# carry out column assays on array s
      prob<- ifelse(sum(array.yk[,j])>0, Se[1], 1-Sp[1]) # compute probability of a positive assay on column j of array s
      g<- rbinom(n=1,size=1,prob=prob)			# generate random assay on column j of array s
      b[j]<-g						# store randomly generated assay on column j of array s in position j of the vector b
      Z[J,1]<-g
      Z[J,2]<-cj
      Z[J,3]<-1
      Z[J,4:(cj+3)]<-array.id[,j]
      Y[array.id[,j],4]<-J
      J<-J+1
    }

    if(sum(a)>0 & sum(b)>0){				# if at least one row and at least one column of array s tested positive...
      array.yk1<-as.matrix(array.yk[(a==1),(b==1)]) 	# pull true statuses of individuals at intersections of rows and columns with positive assays
      array.id1<-as.matrix(array.id[(a==1),(b==1)]) 	# pull individual IDs at intersections of rows and columns with positive assays
      for(i in 1:nrow(array.yk1)){				# carry out individual assays on individuals at intersections of rows and columns with positive assays.
        for(j in 1:ncol(array.yk1)){
          prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])	# get probability of positive assay on individual i,j in intersection
          g<- rbinom(n=1,size=1,prob=prob)			# generate a random assay result for this individual using this probability
          Z[J,1]<-g 						# store this random assay result in row J of Z, column 1
          Z[J,2]<-1 						# store a 1 in row J of Z, column 1, to indicate that this assay was done on 1 individual
          Z[J,3]<-2						# store a 2 in row J of Z, column 2, to be used later to reference Se[2] and Sp[2]
          Z[J,4]<-array.id1[i,j]				# store in row J of Z, column 4, the ID of the individual to whom this assay corresponds
          Y[array.id1[i,j],2]<-3				# store a 3 in the row of Y corresponding to this individual to indicate that this is the third test in which the individual was assayed
          Y[array.id1[i,j],5]<-J				# store J, the index of the assay, in the row of Y corresponding to this individual
          D[array.id1[i,j]] <- g 				#*** store final diagnosis for this individual in the row of D corresponding to this individual
          J<-J+1						# increment assay index J
        }}}

    if(sum(a)>0 & sum(b)==0){				# if no columns but at least one row of array s tested positive...
      array.yk1<-as.matrix(array.yk[(a==1),])		# pull true statuses of individuals in rows of array s which tested positive
      array.id1<-as.matrix(array.id[(a==1),])		# pull IDs of individuals in rows of array s which tested positive
      for(i in 1:nrow(array.yk1)){				# carry out individual assays on individuals in rows of array s which tested positive
        for(j in 1:ncol(array.yk1)){
          prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])# get probability of a positive assay for this individual
          g<- rbinom(n=1,size=1,prob=prob)			# generate random assay result for this individual using this probability
          Z[J,1]<-g 						# store this random assay result in row J of Z, column 1
          Z[J,2]<-1 						# store a 1 in row J of Z, column 2, to indicate that this assay was carried out on 1 individual
          Z[J,3]<-2						# store a 2 in row J of Z, column 3, to be used later to reference Se[2] and Sp[2]
          Z[J,4]<-array.id1[i,j]				# store in row J of Z, column 4, the ID of the individual to whom this assay corresponds
          Y[array.id1[i,j],2]<-3				# store a 3 in the row of Y corresponding to this individual to indicate that this is the third test in which the individual was assayed
          Y[array.id1[i,j],5]<-J				# store J, the index of the assay, in the row of Y corresponding to this individual
          D[array.id1[i,j]] <- g				#*** store final diagnosis for this individual in the row of D corresponding to this individual
          J<-J+1						# increment assay index J
        }}}

    if(sum(a)==0 & sum(b)>0){				# if no rows but at least one column of array s tested positive...
      array.yk1<-as.matrix(array.yk[,(b==1)])
      array.id1<-as.matrix(array.id[,(b==1)])
      for(i in 1:nrow(array.yk1)){
        for(j in 1:ncol(array.yk1)){
          prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
          g<- rbinom(n=1,size=1,prob=prob)
          Z[J,1]<-g
          Z[J,2]<-1
          Z[J,3]<-2
          Z[J,4]<-array.id1[i,j]
          Y[array.id1[i,j],2]<-3
          Y[array.id1[i,j],5]<-J
          D[array.id1[i,j]] <- g				#*** store final diagnosis for this individual in the row of D corresponding to this individual
          J<-J+1
        }}}

  }

  J<-J-1
  Z<-Z[1:J,]

  return(list("Z"=Z, "Y"=Y, "D"=D))
}



#' Approximates the conditional expectations of individual disease statuses with Gibbs sampling.
#'
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param b Parameter values at which to compute the conditional expectations.
#' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
#' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
#' @param GI The length of the Gibbs sampling Markov chain.
#' @return The vector of conditional expectations.
#'
#' This function uses a Gibbs sampler to appriximate the conditional expectation of 
#' each individual's disease status, conditional on the observed assay data and the disease 
#' statuses of all other individuals. This function is used in the EM algorithm
#' performed by the functions \code{mlegt}, \code{enetgt}, \code{enetgt.grid}, and 
#' \code{enetgt.grid.0} under array testing.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' grouplasso_gt_data <- get_grouplasso_gt_data(100)
#' EY <- EYapprox(Z = grouplasso_gt_data$Z,
#'                Y = grouplasso_gt_data$Y,
#'                X = grouplasso_gt_data$X,
#'                b = grouplasso_gt_data$beta,
#'                Se = grouplasso_gt_data$Se,
#'                Sp = grouplasso_gt_data$Sp,
#'                GI = 10000)
#' @export
EYapprox <- function(Z,Y,X,b,Se,Sp,GI=5000){
  n <- dim(Y)[1]
  na <- length(Se)
  # The E-Step
  p <- logit( X %*% b)
  Y[,1] <- rbinom(n,1,p) # this is the initial set of Y values for the Gibbs
  W <- EYgibbs(n,p,Y,Z,Se,Sp,na,GI)
  EY <- W/GI
  return(EY)
}

#' Pulls individual diagnoses from group testing data if available
#'
#' @param \code{Z} output from one of the functions \code{individial.assay.gen()},\code{masterpool.assay.gen()},\code{dorfman.assay.gen()}, or \code{array.assay.gen()}.
#' @param \code{Y} output from one of the functions \code{individial.assay.gen()},\code{masterpool.assay.gen()},\code{dorfman.assay.gen()}, or \code{array.assay.gen()}.
#' @return a vector of 0s and 1s which are the individual diagnoses, NULL if \code{Z} and \code{Y} come from \code{masterpool.assay.gen()}.
#'
#' This function pulls the individual diagnoses from
#' matrices Z and Y when individual diagnoses are available.
#' So for masterpool testing, NULL will be returned.
#' @export
pull.diagnoses <- function(Z,Y)
{
  if( (sum(Z[,2] > 1) == nrow(Z)) & sum(Z[,1]==1) > 0  ) # check to see if masterpool testing has been done.
  {	# all tests on pools of more than one individual? and all tests not negative?

    print("no individual diagnoses available")
    return(NULL)

  }

  # Only individuals who are tested individually can be diagnosed as positive.
  # This is true in Dorfman and Array testing.

  Z.ind.assays <- Z[which(Z[,2]==1),]
  pos.ind <- Z.ind.assays[Z.ind.assays[,1]==1,4]
  D <- numeric(nrow(Y))
  D[pos.ind] <- 1

  return(D)

}


#' Computes conditional expectations of individual disease statuses for individual, master pool, or Dorfman testing
#' 
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param b Parameter values at which to compute the conditional expectations.
#' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
#' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
#' @return The vector of conditional expectations.
#' 
#' This function computes the conditional expectations of each individual
#' disease status, conditional on the observed assay data and the disease
#' statuses of all other individuals.
#' 
#' @examples
#' # generate individual covariate values and disease statuses
#' grouplassogt2pop_data <- get_grouplassogt2pop_data( n1 = 400, n2 = 600)
#' 
#' oldEY <- EYexact_R(Z = grouplassogt2pop_data$Z1,
#'                Y = grouplassogt2pop_data$Y1,
#'                X = grouplassogt2pop_data$X1,
#'                b = rep(1,ncol(grouplassogt2pop_data$X1)),
#'                Se = grouplassogt2pop_data$Se1,
#'                Sp = grouplassogt2pop_data$Sp1)
EYexact_R <- function(Z,Y,X,b,Se,Sp)
{
  # get sample size
  n <- nrow(Y)
  
  # create empty vector in which to store expected values
  EY <- numeric(n)
  
  # get the probs
  p_all <- as.numeric(1/(1 + exp( - X %*% b)))
  
  # get those tested in only one pool
  involved_in_only_one_assay <- c(1:n)[Y[,2] == 1]
  
  while(length(involved_in_only_one_assay) > 0)
  {
    # in which assay was he/she involved?
    assay_in_which_involved <- Y[involved_in_only_one_assay[1],3]
    
    # all involved in this assay
    group_long <- Z[assay_in_which_involved,-c(1,2,3)]
    group <- group_long[group_long != -99]
    
    # get probs for the group
    p_group <- p_all[group]
    
    # get Se and Sp for this assay
    which_SeSp <- Z[assay_in_which_involved,3]
    Sej <- Se[which_SeSp]
    Spj <- Sp[which_SeSp]
    
    # prepare to compute the group probabilities
    prob_whole_group_negative <- prod(1 - p_group)
    zj <- Z[assay_in_which_involved,1]
    Aj <- (Sej*zj + (1-Sej)*(1-zj))
    # Bj <- ((1-Spj)*zj + Spj*(1-zj) - (Sej*zj + (1-Sej)*(1-zj)) ) * prob_whole_group_negative / (1- p_group) + (Sej*zj + (1-Sej)*(1-zj))
    Bj <- ((1-Spj)*zj + Spj*(1-zj) - Aj  ) * prob_whole_group_negative / (1 - p_group) + Aj 
    
    # compute group probabilities
    EY[group] <- Aj * p_group /( Aj * p_group + Bj * (1 - p_group) )
    
    # remove the "processed" individuals from the waiting list
    involved_in_only_one_assay <- involved_in_only_one_assay[ (involved_in_only_one_assay %in% group) == FALSE]
    
  }
  
  # get those involved in two assays
  involved_in_two_assays <- c(1:n)[Y[,2] == 2]
  
  # now for positive Dorfman pools!
  while(length(involved_in_two_assays) > 0)
  {
    
    # in which assays was he/she involved?
    assays_in_which_involved <- Y[involved_in_two_assays[1],3:4]
    
    # all individuals involved in these assays
    group_long <- as.numeric(Z[assays_in_which_involved,-c(1,2,3)])
    group <- unique(group_long[group_long != -99])
    
    # get probs for the group
    p_group <- p_all[group]
    
    # get which assay is the group assay
    group_assay <- Y[group[1],3]
    
    # get individual assays
    individual_assays <- Y[group,4]
    
    # get Se and Sp for the group assay
    which_SeSp <- Z[group_assay,3]
    Sej <- Se[which_SeSp]
    Spj <- Sp[which_SeSp]
    
    # get Se and Sp for the individual assays
    which_SeSp <- Z[individual_assays,3]
    Sej_i <- Se[which_SeSp]
    Spj_i <- Sp[which_SeSp]
    
    # retesting assays
    U <- Z[individual_assays,1]
    
    # prepare mig matrices over which to take sums
    cj <- length(group)
    YY <- expand.grid(rep(list(0:1),cj))
    
    # prepare to compute the group probabilities
    Aj <- numeric(cj)
    Bj <- numeric(cj)
    
    for( k in 1:(2^cj)){
      
      aa <- (Sej_i*U + ( 1 - Sej_i)*(1-U))^YY[k,]
      bb <- ((1-Spj_i)*U + Spj_i*(1-U))^(1 - YY[k,])
      cc <- p_group^YY[k,] * (1 - p_group)^(1-YY[k,])
      
      to_multiply <- as.numeric( aa * bb * cc )
      
      for(i in 1:cj){
        
        prd <- prod(to_multiply[-i])
        max_Y <- max(YY[k,-i])
        Aj[i] <- Aj[i] + prd
        Bj[i] <- Bj[i] + (Sej*max_Y + (1-Spj)*(1 - max_Y)) * prd
        
      }
      
    }
    
    Aj <- Sej * (Sej_i*U + ( 1 - Sej_i)*(1-U)) * Aj
    Bj <- ((1-Spj_i)*U + Spj_i*(1-U)) * Bj
    
    # compute group probabilities
    EY[group] <- Aj * p_group /( Aj * p_group + Bj * (1 - p_group) )
    
    involved_in_two_assays <- involved_in_two_assays[ (involved_in_two_assays %in% group) == FALSE]
    
  }
  
  return(EY)
  
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
                 8,10,
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


#' Generate two data sets for semiparametric additive modeling with group testing responses and some common covariates
#'
#' @param n1 the sample size for the first data set
#' @param n2 the sample size for the second data set
#' @param response a character string indicating the type of response.  Can be \code{"continuous"}, \code{"binary"}, or \code{"gt"}
#' @return a list containing the data
#' @export
get_semipadd2pop_data <- function(n1,n2,response,model = 1,int = FALSE){
  
  p1 <- 6
  q1 <- 4
  zeta1 <- 3/20
  zeta2 <- 10/20
  W1 <- cbind(corrBern(n1,probs=c(1:p1)/(2*p1),Rho = zeta1^abs(outer(1:p1,1:p1,"-"))))
  X1 <- (corrUnif(n1,Rho = zeta2^abs(outer(1:q1,1:q1,"-")))-.5)*5
  
  XX1 <- cbind(1,W1[,c(1,2)],X1[,c(1,2)],W1[,-c(1,2)],X1[,-c(1,2)])
  nonparm1 <- c(0,rep(0,2),rep(1,2),rep(0,p1 - 2),rep(1,q1 - 2))
  pp1 <- ncol(XX1)
  
  # set up the true functions
  f1 <- vector("list",pp1)
  f1[[1]] <- function(x){0} # intercept
  f1[[2]] <- function(x){x*2}
  f1[[3]] <- function(x){x*0}
  f1[[4]] <- function(x){-2*sin(x*2)}
  
  if(model == 1){
    
    f1[[5]] <- function(x){ x }
    
  } else if(model == 2){
    
    f1[[5]] <- function(x){ - x }
    
  }
  
  f1[[6]] <- function(x){x*2}
  f1[[7]] <- function(x){x*0}
  f1[[8]] <- function(x){0*x}
  f1[[9]] <- function(x){0*x}
  f1[[10]] <- function(x){0*x}
  f1[[11]] <- function(x){(exp(-x)-2/5*sinh(5/2))/2}
  
  # record coefficients for covariates to be fit parametrically
  beta1 <- c(-5,2,0,NA,NA,2,0,0,0,NA,NA)
  
  f1.design <- matrix(0,n1,pp1)
  for(j in 1:pp1){
    if(nonparm1[j]==1){
      f1.design[,j] <- f1[[j]](XX1[,j]) - mean(f1[[j]](XX1[,j]))
    } else {
      f1.design[,j] <- f1[[j]](XX1[,j])
    }
  }
  
  # generate second data set
  p2 <- 3
  q2 <- 5
  zeta1 <- 3/20
  zeta2 <- 10/20
  W2 <- cbind(corrBern(n2,probs=c(1:p2)/(2*p2),Rho = zeta1^abs(outer(1:p2,1:p2,"-"))))
  X2 <- (corrUnif(n2,Rho = zeta2^abs(outer(1:q2,1:q2,"-")))-.5)*5
  X2[,1] <- X2[,1] - 1 # impose different supports for some covariates
  X2[,2] <- X2[,2] + 2
  
  XX2 <- cbind(1,W2[,c(1,2)],X2[,c(1,2)],W2[,-c(1,2)],X2[,-c(1,2)])
  nonparm2 <- c(0,rep(0,2),rep(1,2),rep(0,p2 - 2),rep(1,q2 - 2))
  pp2 <- ncol(XX2)
  
  # set up the true functions
  f2 <- vector("list",pp2)
  f2[[1]] <- function(x){ -1 } # intercept
  f2[[2]] <- function(x){ 2 * x  }
  f2[[3]] <- function(x){ 0 * x}
  f2[[4]] <- function(x){ -2*sin(x*2) }
  f2[[5]] <- function(x){ x }
  f2[[6]] <- function(x){ -1 * x }
  f2[[7]] <- function(x){ x^2 - 25/12 }
  f2[[8]] <- function(x){ 0 * x }
  f2[[9]] <- function(x){ 0 * x }
  
  # record coefficients for covariates to be fit parametrically
  beta2 <- c(-4,2,0,NA,NA,-1,NA,NA,NA)
  
  f2.design <- matrix(0,n2,pp2)
  for(j in 1:pp2){
    
    if(nonparm2[j]==1){
      
      f2.design[,j] <- f2[[j]](XX2[,j]) - mean(f2[[j]](XX2[,j]))
      
    } else {
      
      f2.design[,j] <- f2[[j]](XX2[,j])
      
    }
  }
  
  
  
  if( int == TRUE){
    
    
    int1 <- matrix(c(4,2,
                     3,2,
                     3,6),ncol=2,byrow=TRUE)
    
    int2 <- matrix(c(4,2,
                     2,3,
                     6,7),ncol=2,byrow=TRUE)
    
    for( k in 1:nrow(int1)){
      
      if(sum(nonparm1[int1[k,]]) == 1 ){ # int between parametric and nonparametric like f(x)*w
        
        which_nonparm <- int1[k,][which(nonparm1[int1[k,]] == 1)]
        which_parm <- int1[k,][which(nonparm1[int1[k,]] == 0)] 
        
        int.effect1.uncent <- - XX1[,which_nonparm] * XX1[,which_parm] # linear interaction effect modeled nonparametrically
        int.effect1 <- int.effect1.uncent - mean(int.effect1.uncent)
        
      } else if(sum(nonparm1[int1[k,]]) == 0){
        
        int.effect1 <- XX1[,int1[k,1]]*XX1[,int1[k,2]]
        
        
      }
      
      f1.design <- cbind(f1.design,int.effect1)  
      
    }
    
    
    for( k in 1:nrow(int2)){
      
      if(sum(nonparm2[int2[k,]]) == 1 ){ # int between parametric and nonparametric like f(x)*w
        
        which_nonparm <- int2[k,][which(nonparm2[int2[k,]] == 1)]
        which_parm <- int2[k,][which(nonparm2[int2[k,]] == 0)] 
        
        int.effect2.uncent <- - XX2[,which_nonparm] * XX2[,which_parm] # linear interaction effect modeled nonparametrically
        int.effect2 <- int.effect2.uncent - mean(int.effect2.uncent)
        
      } else if(sum(nonparm2[int2[k,]]) == 0){
        
        int.effect2 <- XX2[,int2[k,1]]*XX2[,int2[k,2]]
        
      }
      
      f2.design <- cbind(f2.design,int.effect2)  
      
    }
    
    
    
  } else {
    
    int1 <- NULL
    int2 <- NULL
    
  }
  
  nCom <- 4
  
  if( response == "continuous")
  {
    
    Y1 <- apply(f1.design,1,sum) + rnorm(n1,0,1)
    Y2 <- apply(f2.design,1,sum) + rnorm(n2,0,1)
    
  } else if(response == "binary"){
    
    P1 <- logit(apply(f1.design,1,sum))
    Y1 <- rbinom(n1,1,P1)
    
    P2 <- logit(apply(f2.design,1,sum))
    Y2 <- rbinom(n2,1,P2)
    
  } else if( response == "gt"){
    
    # generate true response values
    P1 <- logit(apply(f1.design,1,sum))
    Y1.true <- rbinom(n1,1,P1)
    
    # generate dorfman testing outcomes
    Se1 <- c(.98,.96)
    Sp1 <- c(.97,.99)
    assay1.out <- dorfman.assay.gen(Y1.true,Se1,Sp1,cj=4)
    Y1 <- list( A = assay1.out$Z,
                I = assay1.out$Y,
                Se = Se1,
                Sp = Sp1,
                E.approx = FALSE)
    
    # generate true response values
    P2 <- logit(apply(f2.design,1,sum))
    Y2.true <- rbinom(n2,1,P2)
    
    # generate individual testing outcomes
    Se2 <- .96
    Sp2 <- .99
    assay2.out <- individual.assay.gen(Y2.true,Se2,Sp2,cj=1)
    Y2 <- list( A = assay2.out$Z,
                I = assay2.out$Y,
                Se = Se2,
                Sp = Sp2,
                E.approx = FALSE)
    
    
  } else {
    
    stop("invalid response type")
    
  }
  
  output <- list(X1 = XX1,
                 nonparm1 = nonparm1,
                 f1 = f1,
                 X2 = XX2,
                 nonparm2 = nonparm2,
                 f2 = f2,
                 Y1 = Y1,
                 Y2 = Y2,
                 nCom = nCom,
                 beta1 = beta1,
                 beta2 = beta2,
                 model = model,
                 int1 = int1,
                 int2 = int2)
  
  return(output)
  
}



#' Generate two data sets for semiparametric additive modeling with continuous responses and some common covariates
#'
#' @param n1 the sample size for the first data set
#' @param n2 the sample size for the second data set
#' @return a list containing the data
#' @export
get_semipadd2pop_simulation_data <- function(n1,n2,some.opposite = TRUE, binary.response = FALSE){

  
  # generate data set 1
  p1 <- 100
  Rho1 <- (.5)^abs(outer(1:p1,1:p1,"-"))
  W1 <- matrix(rnorm(n1*p1),n1,p1) %*% chol(Rho1)
  q1 <- 4
  X1 <- (corrUnif(n1,Rho = (.5)^abs(outer(1:q1,1:q1,"-")))-.5)*5
  
  XX1 <- cbind(1,W1,X1)
  nonparm1 <- c(0,rep(0,p1),rep(1,q1))

  # set up the true functions
  f1 <- vector("list",105)
  f1[[102]] <- function(x){-2*sin(x*2)}
  f1[[103]] <- function(x){x}
  f1[[104]] <- function(x){0*x}
  f1[[105]] <- function(x){(exp(-x)-2/5*sinh(5/2))/2}
  
  # record coefficients for covariates to be fit parametrically
  beta1 <- c(rep(1,5),rep(1,5),rep(0,90))
  
  # generate responses
  if(binary.response){
    
    prob1 <- logit(f1[[102]](XX1[,102]) + f1[[103]](XX1[,103]) + f1[[104]](XX1[,104]) + f1[[105]](XX1[,105]) + W1 %*% beta1)
    
    Y1 <- rbinom(n1,1,prob = prob1)
    
  } else {
    
    Y1 <- f1[[102]](XX1[,102]) + f1[[103]](XX1[,103]) + f1[[104]](XX1[,104]) + f1[[105]](XX1[,105]) + W1 %*% beta1 + rnorm(n1)
  
  }
  
  # generate data set 2
  p2 <- 100
  Rho2 <- (.5)^abs(outer(1:p2,1:p2,"-"))
  W2 <- matrix(rnorm(n2*p2),n2,p2) %*% chol(Rho2)
  q2 <- 5
  X2 <- (corrUnif(n2,Rho = (.5)^abs(outer(1:q2,1:q2,"-")))-.5)*5
  X2[,1] <- X2[,1] - 1 # impose different supports for some covariates
  X2[,2] <- X2[,2] + 2
  
  XX2 <- cbind(1,W2,X2)
  nonparm2 <- c(0,rep(0,p2),rep(1,q2))
  
  # set up the true functions
  f2 <- vector("list",106)
  f2[[102]] <- function(x){-2*sin(x*2)}
  if(some.opposite){
    
    f2[[103]] <- function(x){-x}
    
  } else {
    
    f2[[103]] <- function(x){x}
    
  }
  f2[[104]] <- function(x){x^2 - 25/12}
  f2[[105]] <- function(x){0*x}
  f2[[106]] <- function(x){0*x}
  
  # record coefficients for covariates to be fit parametrically
  if(some.opposite){
    
    beta2 <- c(rep(1,5),rep(-1,5),rep(0,90))
     
  } else {
    
    beta2 <- c(rep(1,5),rep(1,5),rep(0,90))
    
  }
  
  # generate responses
  if(binary.response){
    
    prob2 <- logit(f2[[102]](XX2[,102]) + f2[[103]](XX2[,103]) + f2[[104]](XX2[,104]) + W2 %*% beta2)
    
    Y2 <- rbinom(n2,1,prob = prob2)
    
  } else {
    
    Y2 <- f2[[102]](XX2[,102]) + f2[[103]](XX2[,103]) + f2[[104]](XX2[,104]) + W2 %*% beta2 + rnorm(n2)
    
  }
  
  nCom <- 102
  
  output <- list(X1 = XX1,
                 nonparm1 = nonparm1,
                 f1 = f1,
                 X2 = XX2,
                 nonparm2 = nonparm2,
                 f2 = f2,
                 Y1 = Y1,
                 Y2 = Y2,
                 nCom = nCom,
                 beta1 = c(NA,beta1,rep(NA,q1)),
                 beta2 = c(NA,beta2,rep(NA,q2)),
                 some.opposite = some.opposite,
                 binary.response = binary.response)
  
  return(output)
  
}


