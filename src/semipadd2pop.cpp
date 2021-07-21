# include <RcppArmadillo.h>
// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' The soft-thresholding function for a scalar input
//'
//' @param z the value to which to apply soft-thresholding
//' @param a the threshold
//' @return the value of the soft-thresholding function
//'
//' @examples
//' z <- 3
//' a <- 2
//' SoftThresh(z,a)
// [[Rcpp::export]]
float SoftThresh_scalar(float z, float a){
  
  float sthz;
  
  if( z > a)
  {
    
    sthz = z - a;
    
  } else if( z < - a){
    
    sthz = z + a;
    
  } else {
    
    sthz = 0;
    
  }
  
  return sthz;
  
}

//' Minimize l2-penalized quadratic function
//'
//' @param h a vector
//' @param L a matrix with number of rows equal to the length of h
//' @param lambda a value greater than zero giving the strength of the penalty
//' @param evals the eigenvalues of \eqn{L^TL}
//' @param evecs the eigenvectors of \eqn{L^TL}
//' @return Returns the unique minimizer of \deqn{(1/2) \|h - L \beta\|_2^2  + \lambda * \|\beta\|_2}
//'
//' See Theorem 2 of Foygel, Rina, and Mathias Drton. "Exact block-wise optimization in group lasso and sparse group lasso for linear regression." arXiv preprint arXiv:1010.3320 (2010).
//'
//' @examples
//' # generate an h and L
//' h <- rnorm(100)
//' L <- matrix(rnorm(100*10),100,10)
//' lambda <- 1
//'
//' # get eigendecomposition of t(L) %*% L
//' LtL <- t(L) %*% L
//' eigen.out <- eigen(LtL)
//' evals <- eigen.out$values
//' evecs <- t(eigen.out$vectors)
//'
//' # find minimizer
//' FoygelDrton_Armadillo(h,L,lambda,evals,evecs)
//'
//' # compare to using optim() to minimize the same function
//' obj <- function(beta,L,h,lambda){
//'  val <- (1/2) * sum(  (h - L %*% beta )^2 ) + lambda * sqrt( sum(beta^2))
//'  return(val)
//' }
//' optim(par=rep(0,d),obj,L = L, h = h, lambda = lambda)$par
// [[Rcpp::export]]
arma::colvec FoygelDrton_Armadillo(arma::colvec h, arma::mat L, double lambda, arma::colvec evals, arma::mat evecs) {
  
  arma::colvec v = evecs * arma::trans(L) * h;
  
  double r0, conv, f0, df0, r;
  
  r = 0.1;
  conv = 1;
  
  while(conv > 1e-5){
    
    r0 = r;
    f0 = accu(pow(v,2) / pow(evals * r0 + lambda,2)) - 1;
    df0 = - 2 * accu(  evals % pow(v,2) / pow(evals * r0 + lambda, 3));
    r = std::max(r0 - f0 / df0 , 1e-5);
    conv = std::fabs(r - r0);
    
  }
  
  arma::colvec beta = evecs.t() * arma::diagmat( 1 / (evals + lambda / r)) * v;
  
  return beta;
  
}


//' Minimize the objective function of the group lasso problem with a continuous response
//'
//' @param Y the response vector
//' @param X matrix containing the design matrices
//' @param groups a vector of integers indicating to which group each covariate belongs
//' @param lambda the level of sparsity penalization
//' @param w vector of group-specific weights for different penalization of groups
//' @param eigen a list of eigen info on groups
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta_init optional starting value for beta
//' @return Returns the minimizer of the group lasso objective function
//'
//' @examples
//' data <- get_grouplasso_data(n = 500,response = "continuous")
//' 
//' grouplasso_linreg.out <- grouplasso_linreg(rY = data$Y,
//'                                            rX = data$X,
//'                                            groups = data$groups,
//'                                            lambda = 10,
//'                                            w = data$w,
//'                                            tol = 1e-4,
//'                                            maxiter = 500)
// [[Rcpp::export]]
List grouplasso_linreg_slower(NumericVector rY,
                              NumericMatrix rX,
                              IntegerVector groups,
                              float lambda,
                              NumericVector w,
                              float tol,
                              int maxiter,
                              NumericVector beta_init = NumericVector::create()){
                       
  
  // get dimensions
  int n = rY.size();
  int q = rX.ncol();
  
  int i,j,k,l,i2,k2;
  
  int g = max(groups);
  
  LogicalVector got_eigen(g);
  
  IntegerVector d(g);
  
  IntegerVector is_grp_singleton(g);
  IntegerVector grp_begin_ind(g);
  
  
  // identify singleton and non-singleton groups
  j = 1;
  for(i = 0; i < q; i++){
    
    if(groups[i] != j){
      
      if( d[j-1] == 1){
        
        is_grp_singleton[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d[j-1]++;
    
  }
  
  grp_begin_ind[0] = 0;
  for(j = 1; j < g;j++){
    
    grp_begin_ind[j] = sum(d[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y(rY.begin(),n,false);
  arma::colvec Yj(n);
  arma::colvec rj(n);
  arma::mat X(rX.begin(),n,q,false);
  arma::colvec b = arma::zeros(q);
  
  // Set up initial values if provided
  if( beta_init.size() == q ){
    
    arma::colvec b_silly_copy(beta_init.begin(),q,false);
    b = b_silly_copy;
    
  }
  
  arma::colvec b_0 = b;
  
  arma::colvec Zj_tilde(n);
  arma::colvec h(n);
  arma::mat LtL1(max(d),max(d));
  
  
  arma::field<arma::vec> eigval(g);
  arma::field<arma::mat> eigvec(g);
  arma::field<arma::mat> cholLtL(g);
  
  // define other floats
  float sxj, xjrj_scalar, xjrj_norm;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b_0 = b;
    
    // go through groups
    for( j = 0; j < g ; j++){
      
      // first and last columns of X belonging to group
      i = grp_begin_ind[j];
      k = i + d[j] - 1;
    
      // get predicted values of Y without effects of current covariate
      Yj = arma::zeros(n);
      for( l = 0 ; l < g ; l++)
      {
        
        if(l == j) continue;
        
        i2 = grp_begin_ind[l];
        k2 = i2 + d[l] - 1;
        
        Yj = Yj + X.cols(i2,k2) * b.rows(i2,k2);
        
      }
    
      // get the residuals
      rj = Y - Yj;
      
      if(d[j] == 1){
        
        sxj = arma::as_scalar(X.cols(i,i).t() * X.cols(i,i));
        xjrj_scalar = arma::as_scalar(X.cols(i,i).t() * rj);
        b(i) = SoftThresh_scalar(xjrj_scalar,lambda * w[j] / 2) / sxj;
      
      } else {
        
        xjrj_norm = sqrt(arma::accu(pow(trans(X.cols(i,k)) * rj,2)));
        
        if( xjrj_norm < lambda * w[j]/2){
          
          b.rows(i,k) = arma::zeros(k-i+1);
          
        } else {
          
          if(got_eigen[j] == false){
            
            arma::eig_sym(eigval(j),eigvec(j),X.cols(i,k).t() * X.cols(i,k));
            got_eigen[j] = true;
            
          } 
          
          b.rows(i,k) = FoygelDrton_Armadillo(rj, X.cols(i,k), lambda * w[j], eigval(j), eigvec(j).t());
          
        }
        
      }
      
    }
    
    if(any(abs(b - b_0) > tol)){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta.hat") = b,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}


//' Minimize the objective function of the group lasso problem with a continuous response
//'
//' @param Y the response vector
//' @param X matrix containing the design matrices
//' @param groups a vector of integers indicating to which group each covariate belongs
//' @param lambda the level of sparsity penalization
//' @param w vector of group-specific weights for different penalization of groups
//' @param eigen a list of eigen info on groups
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta_init optional starting value for beta
//' @return Returns the minimizer of the group lasso objective function
//'
//' @examples
//' data <- get_grouplasso_data(n = 500,response = "continuous")
//' 
//' grouplasso_linreg.out <- grouplasso_linreg(rY = data$Y,
//'                                            rX = data$X,
//'                                            groups = data$groups,
//'                                            lambda = 10,
//'                                            w = data$w,
//'                                            tol = 1e-4,
//'                                            maxiter = 500)
// [[Rcpp::export]]
List grouplasso_linreg(NumericVector rY,
                       NumericMatrix rX,
                       IntegerVector groups,
                       double lambda,
                       NumericVector w,
                       double tol,
                       int maxiter,
                       NumericVector beta_init = NumericVector::create()){
                              
  // get dimensions
  int n = rY.size();
  int q = rX.ncol();
  int i,j,k;
  
  int g = max(groups);
  
  LogicalVector got_eigen(g);
  
  IntegerVector d(g);
  
  IntegerVector is_grp_singleton(g);
  IntegerVector grp_begin_ind(g);
  
  // identify singleton and non-singleton groups
  j = 1;
  for(i = 0; i < q; i++){
    
    if(groups[i] != j){
      
      if( d[j-1] == 1){
        
        is_grp_singleton[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d[j-1]++;
    
  }
  
  grp_begin_ind[0] = 0;
  for(j = 1; j < g;j++){
    
    grp_begin_ind[j] = sum(d[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y(rY.begin(),n,false);
  arma::colvec rj(n);
  arma::mat X(rX.begin(),n,q,false);
  arma::colvec b = arma::zeros(q);
  
  // Set up initial values if provided
  if( beta_init.size() == q ){
    
    arma::colvec b_silly_copy(beta_init.begin(),q,false);
    b = b_silly_copy;
    
  }
  
  arma::colvec Y_hat = X * b;
  arma::colvec b_00 = b;
  arma::colvec b_0 = b;
  
  arma::colvec Zj_tilde(n);
  arma::colvec h(n);
  arma::mat LtL1(max(d),max(d));
  
  
  arma::field<arma::vec> eigval(g);
  arma::field<arma::mat> eigvec(g);
  arma::field<arma::mat> cholLtL(g);
  
  // define other floats
  float sxj, xjrj_scalar, xjrj_norm;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b_00 = b;
    
    // go through groups
    for( j = 0; j < g ; j++){
      
      // first and last columns of X belonging to group
      i = grp_begin_ind[j];
      k = i + d[j] - 1;
      
      if(d[j] == 1){
        
        b_0(i) = b(i);
        rj = Y - Y_hat + b(i) * X.cols(i,i);
        
        sxj = arma::as_scalar(X.cols(i,i).t() * X.cols(i,i));
        xjrj_scalar = arma::as_scalar(X.cols(i,i).t() * rj);
        b(i) = SoftThresh_scalar(xjrj_scalar,lambda * w[j] / 2) / sxj;
        
        Y_hat = Y_hat - X.cols(i,i) * ( b_0(i) - b(i) );
        
      } else {
        
        b_0.rows(i,k) = b.rows(i,k);
        rj = Y - Y_hat + X.cols(i,k) * b.rows(i,k);
        
        xjrj_norm = sqrt(arma::accu(pow(trans(X.cols(i,k)) * rj,2)));
        
        if( xjrj_norm < lambda * w[j]/2){
          
          b.rows(i,k) = arma::zeros(k-i+1);
          
        } else {
          
          if(got_eigen[j] == false){
            
            arma::eig_sym(eigval(j),eigvec(j),X.cols(i,k).t() * X.cols(i,k));
            got_eigen[j] = true;
            
          } 
          
          b.rows(i,k) = FoygelDrton_Armadillo(rj, X.cols(i,k), lambda * w[j], eigval(j), eigvec(j).t());
          
        }
        
        Y_hat = Y_hat - X.cols(i,k) * (b_0.rows(i,k) - b.rows(i,k)) ;
        
      }
      
    }
    
    if(any(abs(b - b_00) > tol)){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta.hat") = b,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}


//' Minimize the objective function of the 2-population group lasso problem with a continuous response
//'
//' @param Y1 the continuous response vector of data set 1
//' @param X1 matrix containing the design matrices for data set 1
//' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param Y2 the continuous response vector of data set 2
//' @param X2 matrix containing the design matrices for data set 2
//' @param groups2 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param rho1 weight placed on the first data set
//' @param rho2 weight placed on the second data set
//' @param lambda the level of sparsity penalization
//' @param eta the level of penalization towards model similarity
//' @param w1 group-specific weights for different penalization across groups in data set 1
//' @param w2 group-specific weights for different penalization across groups in data set 2
//' @param w group-specific weights for different penalization toward similarity for different groups
//' @param AA1 a list of the matrices A1j
//' @param AA1 a list of the matrices A2j
//' @param eigen1 a list of eigen info on groups from data set 1
//' @param eigen2 a list of eigen info on groups from data set 2
//' @param Com the indices of the covariate groups which are common in the two datasets
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta1_init optional starting value for beta1
//' @param beta2_init optional starting value for beta2
//' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
//'
//' @examples
//' data <- get_grouplasso2pop_data(n1 = 400, n2 = 600, response = "continuous")
//'   
//' grouplasso2pop_linreg.out <- grouplasso2pop_linreg(rY1 = data$Y1,
//'                                                    rX1 = data$X1,
//'                                                    groups1 = data$groups1,
//'                                                    rY2 = data$Y2,
//'                                                    rX2 = data$X2,
//'                                                    groups2 = data$groups2,
//'                                                    rho1 = 2,
//'                                                    rho2 = 1,
//'                                                    lambda = 1,
//'                                                    eta = 1,
//'                                                    w1 = data$w1,
//'                                                    w2 = data$w2,
//'                                                    w = data$w,
//'                                                    rAA1 = data$AA1,
//'                                                    rAA2 = data$AA2,
//'                                                    rCom = data$Com,
//'                                                    tol = 1e-4,
//'                                                    maxiter = 500)
// [[Rcpp::export]]
List grouplasso2pop_linreg_slower(NumericVector rY1,
                                  NumericMatrix rX1,
                                  IntegerVector groups1,
                                  NumericVector rY2,
                                  NumericMatrix rX2,
                                  IntegerVector groups2,
                                  float rho1,
                                  float rho2,
                                  float lambda,
                                  float eta,
                                  NumericVector w1,
                                  NumericVector w2,
                                  NumericVector w,
                                  List rAA1,
                                  List rAA2,
                                  IntegerVector rCom,
                                  float tol,
                                  int maxiter,
                                  NumericVector beta1_init = NumericVector::create(),
                                  NumericVector beta2_init = NumericVector::create()){
  
  // get dimensions
  int n1 = rY1.size(), n2 = rY2.size();
  int q1 = rX1.ncol(), q2 = rX2.ncol();
  
  int i,j,k,l,m,i2,k2;
  
  int g1 = max(groups1);
  int g2 = max(groups2);
  
  LogicalVector got_eigen1(g1);
  LogicalVector got_eigen2(g2);
  
  IntegerVector d1(g1);
  IntegerVector d2(g2);
  IntegerVector is_grp_singleton1(g1);
  IntegerVector is_grp_singleton2(g2);
  IntegerVector grp_begin_ind1(g1);
  IntegerVector grp_begin_ind2(g2);
  
  // make vector indicating what groups are common
  IntegerVector is_com(std::max(g1,g2));
  
  for(j = 0; j < rCom.size(); j++){
    
    is_com[rCom[j]-1] = 1;
    
  }
  
  // identify singleton and non-singleton groups for data set 1 
  j = 1;
  for(i = 0; i < q1; i++){
    
    if(groups1[i] != j){
      
      if( d1[j-1] == 1){
        
        is_grp_singleton1[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d1[j-1]++;
    
  }
  
  grp_begin_ind1[0] = 0;
  for(j = 1; j < g1;j++){
    
    grp_begin_ind1[j] = sum(d1[Range(0,j-1)]);
    
  }
  
  // identify singleton and non-singleton groups for data set 2
  j = 1;
  for(i = 0; i < q2; i++){
    
    if(groups2[i] != j){
      
      if( d2[j-1] == 1){
        
        is_grp_singleton2[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d2[j-1] ++;
    
  }
  
  grp_begin_ind2[0] = 0;
  for(j = 1; j < g2;j++){
    
    grp_begin_ind2[j] = sum(d2[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y1(rY1.begin(),n1,false);
  arma::colvec Y2(rY2.begin(),n2,false);
  
  arma::colvec Y1j(n1);
  arma::colvec r1j(n1);
  
  arma::colvec Y2j(n2);
  arma::colvec r2j(n2);
  
  arma::mat X1(rX1.begin(),n1,q1,false);
  arma::mat X2(rX2.begin(),n2,q2,false);
  
  arma::colvec b1 = arma::zeros(q1);
  arma::colvec b2 = arma::zeros(q2);
  
  // Set up initial values if provided
  if( beta1_init.size() == q1 ){
    
    arma::colvec b1_silly_copy(beta1_init.begin(),q1,false);
    b1 = b1_silly_copy;
    
  }
  
  if( beta2_init.size() == q2 ){
    
    arma::colvec b2_silly_copy(beta2_init.begin(),q2,false);
    b2 = b2_silly_copy;
    
  }
  
  arma::colvec b1_0 = b1;
  arma::colvec b2_0 = b2;
  
  arma::colvec h1(n1);
  arma::colvec h2(n2);
  
  arma::colvec x1r1j_AAb2(max(d1));
  arma::colvec x2r2j_AAb1(max(d2));
  arma::mat LtL1(max(d1),max(d1));
  arma::mat LtL2(max(d2),max(d2));
  
  arma::field<arma::vec> eigval1(g1);
  arma::field<arma::mat> eigvec1(g1);
  arma::field<arma::mat> cholLtL1(g1);
  
  arma::field<arma::vec> eigval2(g2);
  arma::field<arma::mat> eigvec2(g2);
  arma::field<arma::mat> cholLtL2(g2);
  
  k = max(rCom);
  arma::field<arma::mat> AA1(k);
  arma::field<arma::mat> AA2(k);
  
  for( j = 0; j < k; j++){
    if(is_com[j] == 1){
      
      arma::mat AA1_silly_copy(as<NumericMatrix>(rAA1[j]).begin(),as<NumericMatrix>(rAA1[j]).nrow(),d1[j]);
      AA1(j) = AA1_silly_copy;
      
      arma::mat AA2_silly_copy(as<NumericMatrix>(rAA2[j]).begin(),as<NumericMatrix>(rAA2[j]).nrow(),d2[j]);
      AA2(j) = AA2_silly_copy;
      
    }
    
  }
  
  // define other floats
  float sx1j, x1r1j_scalar, x1r1j_norm;
  float sx1j_AA, x1r1j_AAb2_scalar, x1r1j_AAb2_norm;
  
  float sx2j, x2r2j_scalar, x2r2j_norm;
  float sx2j_AA, x2r2j_AAb1_scalar, x2r2j_AAb1_norm;
  
  float lambda1 = lambda * n1 / rho1;
  float eta1 = eta * n1 / rho1;
  
  float lambda2 = lambda * n2 / rho1;
  float eta2 = eta * n2 / rho2;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b1_0 = b1;
    b2_0 = b2;
    
    // go through groups of data set 1
    for( j = 0; j < g1 ; j++){
      
      // first and last columns of X1 belonging to group
      i = grp_begin_ind1[j];
      k = i + d1[j] - 1;
      
      // get predicted values of Y without effects of current covariate
      Y1j = arma::zeros(n1);
      for( l = 0 ; l < g1 ; l++)
      {
        
        if(l == j) continue;
        
        i2 = grp_begin_ind1[l];
        k2 = i2 + d1[l] - 1;
        
        Y1j = Y1j + X1.cols(i2,k2) * b1.rows(i2,k2);
        
      }
      
      // get residuals
      r1j = Y1 - Y1j;
      
      if(is_com[j] == 1){
        
        // first and last columns of X2 belonging to corresponding group
        l = grp_begin_ind2[j];
        m = l + d2[j] - 1;
        
        if(d1[j] == 1){
          
          x1r1j_AAb2_scalar = arma::as_scalar(X1.cols(i,i).t() * r1j + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m));
          sx1j_AA = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i) + eta1 * w[j] * AA1(j).t() * AA1(j));
          b1(i) = SoftThresh_scalar(x1r1j_AAb2_scalar,lambda1*w1[j]) / sx1j_AA;  
          
        } else {
          
          x1r1j_AAb2.rows(0,d1[j]-1) = X1.cols(i,k).t() * r1j + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m);
          x1r1j_AAb2_norm = sqrt(accu(pow(x1r1j_AAb2.rows(0,d1[j]-1),2)));
          
          if(x1r1j_AAb2_norm < lambda1 * w1[j])
          {
            
            b1.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen1[j] == false){
              
              LtL1.submat(0,0,d1[j]-1,d1[j]-1) = X1.cols(i,k).t() * X1.cols(i,k) + eta1 * w[j] * AA1(j).t() * AA1(j);
              arma::eig_sym(eigval1(j),eigvec1(j),LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              cholLtL1(j) = chol(LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              got_eigen1[j] = true;
              
            }
            
            h1 = arma::inv(cholLtL1(j).t()) * x1r1j_AAb2.rows(0,d1[j]-1);
            
            b1.rows(i,k) = FoygelDrton_Armadillo(h1, cholLtL1(j), lambda1 * w1[j], eigval1(j), eigvec1(j).t());
            
          }
          
        }
        
      } else { // not a common covariate
        
        
        if(d1[j] == 1){
          
          sx1j = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i));
          x1r1j_scalar = arma::as_scalar(X1.cols(i,i).t() * r1j);
          b1(i) = SoftThresh_scalar(x1r1j_scalar,lambda1 * w1[j]) / sx1j;
          
        } else {
          
          x1r1j_norm = sqrt(arma::accu(pow(trans(X1.cols(i,k)) * r1j,2)));
          
          if( x1r1j_norm < lambda1 * w1[j]){
            
            b1.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen1[j] == false){
              
              arma::eig_sym(eigval1(j),eigvec1(j),X1.cols(i,k).t() * X1.cols(i,k));
              got_eigen1[j] = true;
              
            } 
            
            b1.rows(i,k) = FoygelDrton_Armadillo(r1j, X1.cols(i,k), lambda1 * w1[j], eigval1(j), eigvec1(j).t());
            
          }
          
        }
        
      }
      
    }
    
    // go through groups of data set 2
    for( j = 0; j < g2 ; j++){
      
      // first and last columns of X2 belonging to group
      i = grp_begin_ind2[j];
      k = i + d2[j] - 1;
      
      // get predicted values of Y without effects of current covariate
      Y2j = arma::zeros(n2);
      for( l = 0 ; l < g2 ; l++)
      {
        
        if(l == j) continue;
        
        i2 = grp_begin_ind2[l];
        k2 = i2 + d2[l] - 1;
        
        Y2j = Y2j + X2.cols(i2,k2) * b2.rows(i2,k2);
        
      }
      
      // get residuals
      r2j = Y2 - Y2j;
      
      if(is_com[j] == 1){
        
        // first and last columns of X1 belonging to corresponding group
        l = grp_begin_ind1[j];
        m = l + d1[j] - 1;
        
        if(d2[j] == 1){
          
          x2r2j_AAb1_scalar = arma::as_scalar(X2.cols(i,i).t() * r2j + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m));
          sx2j_AA = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i) + eta2 * w[j] * AA2(j).t() * AA2(j));
          b2(i) = SoftThresh_scalar(x2r2j_AAb1_scalar,lambda2*w2[j]) / sx2j_AA;  
          
        } else {
          
          x2r2j_AAb1.rows(0,d2[j]-1) = X2.cols(i,k).t() * r2j + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m);
          x2r2j_AAb1_norm = sqrt(accu(pow(x2r2j_AAb1.rows(0,d2[j]-1),2)));
          
          if(x2r2j_AAb1_norm < lambda2 * w2[j])
          {
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              LtL2.submat(0,0,d2[j]-1,d2[j]-1) = X2.cols(i,k).t() * X2.cols(i,k) + eta2 * w[j] * AA2(j).t() * AA2(j);
              arma::eig_sym(eigval2(j),eigvec2(j),LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              cholLtL2(j) = chol(LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              got_eigen2[j] = true;
              
            }
            
            h2 = arma::inv(cholLtL2(j).t()) * x2r2j_AAb1.rows(0,d2[j]-1);
            b2.rows(i,k) = FoygelDrton_Armadillo(h2, cholLtL2(j), lambda2 * w2[j], eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      } else { // not a common covariate
        
        
        if(d2[j] == 1){
          
          sx2j = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i));
          x2r2j_scalar = arma::as_scalar(X2.cols(i,i).t() * r2j);
          b2(i) = SoftThresh_scalar(x2r2j_scalar,lambda2 * w2[j]) / sx2j;
          
        } else {
          
          x2r2j_norm = sqrt(arma::accu(pow(trans(X2.cols(i,k)) * r2j,2)));
          
          if( x2r2j_norm < lambda2 * w2[j]){
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              arma::eig_sym(eigval2(j),eigvec2(j),X2.cols(i,k).t() * X2.cols(i,k));
              got_eigen2[j] = true;
              
            } 
            
            b2.rows(i,k) = FoygelDrton_Armadillo(r2j, X2.cols(i,k), lambda2 * w2[j], eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      }
      
    }
    
    if(any(abs(b1 - b1_0) > tol) | any(abs(b2 - b2_0) > tol) ){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta1.hat") = b1,
                            Named("beta2.hat") = b2,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}

//' Minimize the objective function of the 2-population group lasso problem with a continuous response
//'
//' @param Y1 the continuous response vector of data set 1
//' @param X1 matrix containing the design matrices for data set 1
//' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param Y2 the continuous response vector of data set 2
//' @param X2 matrix containing the design matrices for data set 2
//' @param groups2 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param rho1 weight placed on the first data set
//' @param rho2 weight placed on the second data set
//' @param lambda the level of sparsity penalization
//' @param eta the level of penalization towards model similarity
//' @param w1 group-specific weights for different penalization across groups in data set 1
//' @param w2 group-specific weights for different penalization across groups in data set 2
//' @param w group-specific weights for different penalization toward similarity for different groups
//' @param AA1 a list of the matrices A1j
//' @param AA1 a list of the matrices A2j
//' @param eigen1 a list of eigen info on groups from data set 1
//' @param eigen2 a list of eigen info on groups from data set 2
//' @param Com the indices of the covariate groups which are common in the two datasets
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta1_init optional starting value for beta1
//' @param beta2_init optional starting value for beta2
//' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
//'
//' @examples
//' data <- get_grouplasso2pop_data(n1 = 400, n2 = 600, response = "continuous")
//'   
//' grouplasso2pop_linreg.out <- grouplasso2pop_linreg(rY1 = data$Y1,
//'                                                    rX1 = data$X1,
//'                                                    groups1 = data$groups1,
//'                                                    rY2 = data$Y2,
//'                                                    rX2 = data$X2,
//'                                                    groups2 = data$groups2,
//'                                                    rho1 = 2,
//'                                                    rho2 = 1,
//'                                                    lambda = 1,
//'                                                    eta = 1,
//'                                                    w1 = data$w1,
//'                                                    w2 = data$w2,
//'                                                    w = data$w,
//'                                                    rAA1 = data$AA1,
//'                                                    rAA2 = data$AA2,
//'                                                    rCom = data$Com,
//'                                                    tol = 1e-4,
//'                                                    maxiter = 500)
// [[Rcpp::export]]
List grouplasso2pop_linreg(NumericVector rY1,
                           NumericMatrix rX1,
                           IntegerVector groups1,
                           NumericVector rY2,
                           NumericMatrix rX2,
                           IntegerVector groups2,
                           float rho1,
                           float rho2,
                           float lambda,
                           float eta,
                           NumericVector w1,
                           NumericVector w2,
                           NumericVector w,
                           List rAA1,
                           List rAA2,
                           IntegerVector rCom,
                           float tol,
                           int maxiter,
                           NumericVector beta1_init = NumericVector::create(),
                           NumericVector beta2_init = NumericVector::create()){
  
  // get dimensions
  int n1 = rY1.size(), n2 = rY2.size();
  int q1 = rX1.ncol(), q2 = rX2.ncol();
  
  int i,j,k,l,m;
  
  int g1 = max(groups1);
  int g2 = max(groups2);
  
  LogicalVector got_eigen1(g1);
  LogicalVector got_eigen2(g2);
  
  IntegerVector d1(g1);
  IntegerVector d2(g2);
  IntegerVector is_grp_singleton1(g1);
  IntegerVector is_grp_singleton2(g2);
  IntegerVector grp_begin_ind1(g1);
  IntegerVector grp_begin_ind2(g2);
  
  // make vector indicating what groups are common
  IntegerVector is_com(std::max(g1,g2));
  
  for(j = 0; j < rCom.size(); j++){
    
    is_com[rCom[j]-1] = 1;
    
  }
  
  // identify singleton and non-singleton groups for data set 1 
  j = 1;
  for(i = 0; i < q1; i++){
    
    if(groups1[i] != j){
      
      if( d1[j-1] == 1){
        
        is_grp_singleton1[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d1[j-1]++;
    
  }
  
  grp_begin_ind1[0] = 0;
  for(j = 1; j < g1;j++){
    
    grp_begin_ind1[j] = sum(d1[Range(0,j-1)]);
    
  }
  
  // identify singleton and non-singleton groups for data set 2
  j = 1;
  for(i = 0; i < q2; i++){
    
    if(groups2[i] != j){
      
      if( d2[j-1] == 1){
        
        is_grp_singleton2[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d2[j-1] ++;
    
  }
  
  grp_begin_ind2[0] = 0;
  for(j = 1; j < g2;j++){
    
    grp_begin_ind2[j] = sum(d2[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y1(rY1.begin(),n1,false);
  arma::colvec Y2(rY2.begin(),n2,false);
  
  arma::colvec Y1j(n1);
  arma::colvec r1j(n1);
  
  arma::colvec Y2j(n2);
  arma::colvec r2j(n2);
  
  arma::mat X1(rX1.begin(),n1,q1,false);
  arma::mat X2(rX2.begin(),n2,q2,false);
  
  arma::colvec b1 = arma::zeros(q1);
  arma::colvec b2 = arma::zeros(q2);
  
  // Set up initial values if provided
  if( beta1_init.size() == q1 ){
    
    arma::colvec b1_silly_copy(beta1_init.begin(),q1,false);
    b1 = b1_silly_copy;
    
  }
  
  if( beta2_init.size() == q2 ){
    
    arma::colvec b2_silly_copy(beta2_init.begin(),q2,false);
    b2 = b2_silly_copy;
    
  }
  
  arma::colvec b1_00 = b1;
  arma::colvec b1_0 = b1;
  arma::colvec b2_00 = b2;
  arma::colvec b2_0 = b2;
  
  arma::colvec Y1_hat = X1 * b1;
  arma::colvec Y2_hat = X2 * b2;
  
  arma::colvec h1(n1);
  arma::colvec h2(n2);
  
  arma::colvec x1r1j_AAb2(max(d1));
  arma::colvec x2r2j_AAb1(max(d2));
  arma::mat LtL1(max(d1),max(d1));
  arma::mat LtL2(max(d2),max(d2));
  
  arma::field<arma::vec> eigval1(g1);
  arma::field<arma::mat> eigvec1(g1);
  arma::field<arma::mat> cholLtL1(g1);
  
  arma::field<arma::vec> eigval2(g2);
  arma::field<arma::mat> eigvec2(g2);
  arma::field<arma::mat> cholLtL2(g2);
  
  k = max(rCom);
  arma::field<arma::mat> AA1(k);
  arma::field<arma::mat> AA2(k);
  
  for( j = 0; j < k; j++){
    if(is_com[j] == 1){
      
      arma::mat AA1_silly_copy(as<NumericMatrix>(rAA1[j]).begin(),as<NumericMatrix>(rAA1[j]).nrow(),d1[j]);
      AA1(j) = AA1_silly_copy;
      
      arma::mat AA2_silly_copy(as<NumericMatrix>(rAA2[j]).begin(),as<NumericMatrix>(rAA2[j]).nrow(),d2[j]);
      AA2(j) = AA2_silly_copy;
      
    }
    
  }
  
  // define other floats
  float sx1j, x1r1j_scalar, x1r1j_norm;
  float sx1j_AA, x1r1j_AAb2_scalar, x1r1j_AAb2_norm;
  
  float sx2j, x2r2j_scalar, x2r2j_norm;
  float sx2j_AA, x2r2j_AAb1_scalar, x2r2j_AAb1_norm;
  
  float lambda1 = lambda * n1 / rho1;
  float eta1 = eta * n1 / rho1;
  
  float lambda2 = lambda * n2 / rho1;
  float eta2 = eta * n2 / rho2;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b1_00 = b1;
    b2_00 = b2;
    
    // go through groups of data set 1
    for( j = 0; j < g1 ; j++){
      
      // first and last columns of X1 belonging to group
      i = grp_begin_ind1[j];
      k = i + d1[j] - 1;
      
      // get residuals
      r1j = Y1 - Y1_hat + X1.cols(i,k) * b1.rows(i,k);
      b1_0.rows(i,k) = b1.rows(i,k);
      
      if(is_com[j] == 1){
        
        // first and last columns of X2 belonging to corresponding group
        l = grp_begin_ind2[j];
        m = l + d2[j] - 1;
        
        if(d1[j] == 1){
          
          x1r1j_AAb2_scalar = arma::as_scalar(X1.cols(i,i).t() * r1j + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m));
          sx1j_AA = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i) + eta1 * w[j] * AA1(j).t() * AA1(j));
          b1(i) = SoftThresh_scalar(x1r1j_AAb2_scalar,lambda1*w1[j]) / sx1j_AA;  
          
        } else {
          
          x1r1j_AAb2.rows(0,d1[j]-1) = X1.cols(i,k).t() * r1j + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m);
          x1r1j_AAb2_norm = sqrt(accu(pow(x1r1j_AAb2.rows(0,d1[j]-1),2)));
          
          if(x1r1j_AAb2_norm < lambda1 * w1[j])
          {
            
            b1.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen1[j] == false){
  
              LtL1.submat(0,0,d1[j]-1,d1[j]-1) = X1.cols(i,k).t() * X1.cols(i,k) + eta1 * w[j] * AA1(j).t() * AA1(j);
              arma::eig_sym(eigval1(j),eigvec1(j),LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              cholLtL1(j) = chol(LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              got_eigen1[j] = true;
  
            }
  
            h1 = arma::inv(cholLtL1(j).t()) * x1r1j_AAb2.rows(0,d1[j]-1);
  
            b1.rows(i,k) = FoygelDrton_Armadillo(h1, cholLtL1(j), lambda1 * w1[j], eigval1(j), eigvec1(j).t());
              
          }
          
        }
        
      } else { // not a common covariate
        
        
          if(d1[j] == 1){
            
            sx1j = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i));
            x1r1j_scalar = arma::as_scalar(X1.cols(i,i).t() * r1j);
            b1(i) = SoftThresh_scalar(x1r1j_scalar,lambda1 * w1[j]) / sx1j;
            
          } else {
            
            x1r1j_norm = sqrt(arma::accu(pow(trans(X1.cols(i,k)) * r1j,2)));
            
            if( x1r1j_norm < lambda1 * w1[j]){
              
              b1.rows(i,k) = arma::zeros(k-i+1);
              
            } else {
              
              if(got_eigen1[j] == false){
                
                arma::eig_sym(eigval1(j),eigvec1(j),X1.cols(i,k).t() * X1.cols(i,k));
                got_eigen1[j] = true;
                
              } 
              
              b1.rows(i,k) = FoygelDrton_Armadillo(r1j, X1.cols(i,k), lambda1 * w1[j], eigval1(j), eigvec1(j).t());
              
            }
            
          }
          
        }
        
        Y1_hat = Y1_hat - X1.cols(i,k) * ( b1_0.rows(i,k)  - b1.rows(i,k));
        
      }
      
    // go through groups of data set 2
    for( j = 0; j < g2 ; j++){
      
      // first and last columns of X2 belonging to group
      i = grp_begin_ind2[j];
      k = i + d2[j] - 1;
      
      
      // get residuals
      r2j = Y2 - Y2_hat + X2.cols(i,k) * b2.rows(i,k);
      b2_0.rows(i,k) = b2.rows(i,k);
      
      if(is_com[j] == 1){
        
        // first and last columns of X1 belonging to corresponding group
        l = grp_begin_ind1[j];
        m = l + d1[j] - 1;
        
        if(d2[j] == 1){
          
          x2r2j_AAb1_scalar = arma::as_scalar(X2.cols(i,i).t() * r2j + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m));
          sx2j_AA = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i) + eta2 * w[j] * AA2(j).t() * AA2(j));
          b2(i) = SoftThresh_scalar(x2r2j_AAb1_scalar,lambda2*w2[j]) / sx2j_AA;  
          
        } else {
          
          x2r2j_AAb1.rows(0,d2[j]-1) = X2.cols(i,k).t() * r2j + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m);
          x2r2j_AAb1_norm = sqrt(accu(pow(x2r2j_AAb1.rows(0,d2[j]-1),2)));
          
          if(x2r2j_AAb1_norm < lambda2 * w2[j])
          {
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              LtL2.submat(0,0,d2[j]-1,d2[j]-1) = X2.cols(i,k).t() * X2.cols(i,k) + eta2 * w[j] * AA2(j).t() * AA2(j);
              arma::eig_sym(eigval2(j),eigvec2(j),LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              cholLtL2(j) = chol(LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              got_eigen2[j] = true;
              
            }
            
            h2 = arma::inv(cholLtL2(j).t()) * x2r2j_AAb1.rows(0,d2[j]-1);
            b2.rows(i,k) = FoygelDrton_Armadillo(h2, cholLtL2(j), lambda2 * w2[j], eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      } else { // not a common covariate
        
        
        if(d2[j] == 1){
          
          sx2j = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i));
          x2r2j_scalar = arma::as_scalar(X2.cols(i,i).t() * r2j);
          b2(i) = SoftThresh_scalar(x2r2j_scalar,lambda2 * w2[j]) / sx2j;
          
        } else {
          
          x2r2j_norm = sqrt(arma::accu(pow(trans(X2.cols(i,k)) * r2j,2)));
          
          if( x2r2j_norm < lambda2 * w2[j]){
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              arma::eig_sym(eigval2(j),eigvec2(j),X2.cols(i,k).t() * X2.cols(i,k));
              got_eigen2[j] = true;
              
            } 
            
            b2.rows(i,k) = FoygelDrton_Armadillo(r2j, X2.cols(i,k), lambda2 * w2[j], eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      }
      
      Y2_hat = Y2_hat - X2.cols(i,k) * ( b2_0.rows(i,k)  - b2.rows(i,k));
      
    }
    
    if(any(abs(b1 - b1_00) > tol) | any(abs(b2 - b2_00) > tol) ){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta1.hat") = b1,
                            Named("beta2.hat") = b2,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}




//' Minimize the objective function of the group lasso problem with a binary response
//'
//' @param Y the binary response vector
//' @param X matrix containing the design matrices
//' @param groups a vector of integers indicating to which group each covariate belongs
//' @param lambda the level of sparsity penalization
//' @param w vector of group-specific weights for different penalization of groups
//' @param eigen a list of eigen info on groups
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta_init optional starting value for beta
//' @return Returns the minimizer of the group lasso objective function
//'
//' @examples
//' data <- get_grouplasso_data(n = 500, response = "binary")
//' 
//' grouplasso_logreg.out <- grouplasso_logreg(rY = data$Y,
//'                                            rX = data$X,
//'                                            groups = data$groups,
//'                                            lambda = 10,
//'                                            w = data$w,
//'                                            tol = 1e-4,
//'                                            maxiter = 500)
// [[Rcpp::export]]
List grouplasso_logreg_slower(NumericVector rY,
                              NumericMatrix rX,
                              IntegerVector groups,
                              float lambda,
                              NumericVector w,
                              float tol,
                              int maxiter,
                              NumericVector beta_init = NumericVector::create()){
                       
                      
  // get dimensions
  int n = rY.size();
  int q = rX.ncol();
  
  int i,j,k;
  
  int g = max(groups);
  
  LogicalVector got_eigen(g);
  
  IntegerVector d(g);
  
  IntegerVector is_grp_singleton(g);
  IntegerVector grp_begin_ind(g);
  

  // identify singleton and non-singleton groups
  j = 1;
  for(i = 0; i < q; i++){
    
    if(groups[i] != j){
      
      if( d[j-1] == 1){
        
        is_grp_singleton[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d[j-1]++;
    
  }
  
  grp_begin_ind[0] = 0;
  for(j = 1; j < g;j++){
    
    grp_begin_ind[j] = sum(d[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y(rY.begin(),n,false);
  arma::colvec P(n);
  arma::mat X(rX.begin(),n,q,false);
  arma::colvec b = arma::zeros(q);
  
  // Set up initial values if provided
  if( beta_init.size() == q ){
    
    arma::colvec b_silly_copy(beta_init.begin(),q,false);
    b = b_silly_copy;
    
  }
  
  arma::colvec b_0 = b;

  arma::colvec Zj_tilde(n);
  arma::colvec h(n);
  arma::mat LtL1(max(d),max(d));
  
  
  arma::field<arma::vec> eigval(g);
  arma::field<arma::mat> eigvec(g);
  arma::field<arma::mat> cholLtL(g);
  
  // define other floats
  float sxj, sxzj_scalar, xzj_norm;
  //float sx1j_AA, sx1z1j_AAb2_scalar, x1z1j_AAb2_norm;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b_0 = b;
  
    // go through groups
    for( j = 0; j < g ; j++){
      
      // update probs
      P = 1/(1 + exp( - X * b  ));
      
      if(any(P > 1-1e-10) | any(P < 1e-10)){
        
        Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
        iter = maxiter;
        
      }
      
      // first and last columns of X belonging to group
      i = grp_begin_ind[j];
      k = i + d[j] - 1;

      if(d[j] == 1){
        
        
        Zj_tilde = X.cols(i,i) * b(i) / 2 + 2 * ( Y - P);
        sxj = arma::as_scalar(X.cols(i,i).t() * X.cols(i,i)) / 4;
        sxzj_scalar = arma::as_scalar(X.cols(i,i).t() * Zj_tilde) / 2;
        b(i) = SoftThresh_scalar(sxzj_scalar,lambda * w[j] / 2) / sxj;
        
      } else {
        
        Zj_tilde = X.cols(i,k) * b.rows(i,k) / 2 + 2 * ( Y - P);
        xzj_norm = sqrt(arma::accu(pow(trans(X.cols(i,k)) * Zj_tilde,2))) / 2;
        
        if( xzj_norm < lambda * w[j]/2){
          
          b.rows(i,k) = arma::zeros(k-i+1);
          
        } else {
          
          if(got_eigen[j] == false){
            
            arma::eig_sym(eigval(j),eigvec(j),X.cols(i,k).t() * X.cols(i,k) / 4);
            got_eigen[j] = true;
            
          } 
          
          b.rows(i,k) = FoygelDrton_Armadillo(Zj_tilde, X.cols(i,k) / 2, lambda * w[j]/2, eigval(j), eigvec(j).t());
          
        }
        
      }
      
    }
    
    if(any(abs(b - b_0) > tol)){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta.hat") = b,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}


//' Minimize the objective function of the group lasso problem with a binary response
//'
//' @param Y the binary response vector
//' @param X matrix containing the design matrices
//' @param groups a vector of integers indicating to which group each covariate belongs
//' @param lambda the level of sparsity penalization
//' @param w vector of group-specific weights for different penalization of groups
//' @param eigen a list of eigen info on groups
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta_init optional starting value for beta
//' @return Returns the minimizer of the group lasso objective function
//'
//' @examples
//' data <- get_grouplasso_data(n = 500, response = "binary")
//' 
//' grouplasso_logreg.out <- grouplasso_logreg(rY = data$Y,
//'                                            rX = data$X,
//'                                            groups = data$groups,
//'                                            lambda = 10,
//'                                            w = data$w,
//'                                            tol = 1e-4,
//'                                            maxiter = 500)
// [[Rcpp::export]]
List grouplasso_logreg(NumericVector rY,
                       NumericMatrix rX,
                       IntegerVector groups,
                       float lambda,
                       NumericVector w,
                       float tol,
                       int maxiter,
                       NumericVector beta_init = NumericVector::create()){
                              
  // get dimensions
  int n = rY.size();
  int q = rX.ncol();
  
  int i,j,k;
  
  int g = max(groups);
  
  LogicalVector got_eigen(g);
  
  IntegerVector d(g);
  
  IntegerVector is_grp_singleton(g);
  IntegerVector grp_begin_ind(g);
  
  
  // identify singleton and non-singleton groups
  j = 1;
  for(i = 0; i < q; i++){
    
    if(groups[i] != j){
      
      if( d[j-1] == 1){
        
        is_grp_singleton[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d[j-1]++;
    
  }
  
  grp_begin_ind[0] = 0;
  for(j = 1; j < g;j++){
    
    grp_begin_ind[j] = sum(d[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y(rY.begin(),n,false);
  arma::colvec P(n);
  arma::mat X(rX.begin(),n,q,false);
  arma::colvec b = arma::zeros(q);
  
  // Set up initial values if provided
  if( beta_init.size() == q ){
    
    arma::colvec b_silly_copy(beta_init.begin(),q,false);
    b = b_silly_copy;
    
  }
  
  arma::colvec log_odds = X * b;
  arma::colvec b_00 = b;
  arma::colvec b_0 = b;
  
  arma::colvec Zj_tilde(n);
  arma::colvec h(n);
  arma::mat LtL1(max(d),max(d));
  
  arma::field<arma::vec> eigval(g);
  arma::field<arma::mat> eigvec(g);
  arma::field<arma::mat> cholLtL(g);
  
  // define other floats
  float sxj, sxzj_scalar, xzj_norm;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b_00 = b;
    
    // go through groups
    for( j = 0; j < g ; j++){
      
      // update probs
      P = 1/(1 + exp( - log_odds  ));
      
      if(any(P > 1-1e-10) | any(P < 1e-10)){
        
        Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
        iter = maxiter;
        
      }
      
      // first and last columns of X belonging to group
      i = grp_begin_ind[j];
      k = i + d[j] - 1;
      
      b_0.rows(i,k) = b.rows(i,k);
      
      if(d[j] == 1){
        
        Zj_tilde = X.cols(i,i) * b(i) / 2 + 2 * ( Y - P);
        sxj = arma::as_scalar(X.cols(i,i).t() * X.cols(i,i)) / 4;
        sxzj_scalar = arma::as_scalar(X.cols(i,i).t() * Zj_tilde) / 2;
        b(i) = SoftThresh_scalar(sxzj_scalar,lambda * w[j] / 2) / sxj;
        
        log_odds = log_odds - X.cols(i,i) * (b_0(i) - b(i));
        
      } else {
        
        Zj_tilde = X.cols(i,k) * b.rows(i,k) / 2 + 2 * ( Y - P);
        xzj_norm = sqrt(arma::accu(pow(trans(X.cols(i,k)) * Zj_tilde,2))) / 2;
        
        if( xzj_norm < lambda * w[j]/2){
          
          b.rows(i,k) = arma::zeros(k-i+1);
          
        } else {
          
          if(got_eigen[j] == false){
            
            arma::eig_sym(eigval(j),eigvec(j),X.cols(i,k).t() * X.cols(i,k) / 4);
            got_eigen[j] = true;
            
          } 
          
          b.rows(i,k) = FoygelDrton_Armadillo(Zj_tilde, X.cols(i,k) / 2, lambda * w[j]/2, eigval(j), eigvec(j).t());
          
        }
        
        log_odds = log_odds - X.cols(i,k) * (b_0.rows(i,k) - b.rows(i,k));
        
      }
      
    }
    
    if(any(abs(b - b_00) > tol)){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta.hat") = b,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}



//' Minimize the objective function of the 2-population group lasso problem with a binary response
//'
//' @param Y1 the binary response vector of data set 1
//' @param X1 matrix containing the design matrices for data set 1
//' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param Y2 the binary response vector of data set 2
//' @param X2 matrix containing the design matrices for data set 2
//' @param groups2 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param rho1 weight placed on the first data set
//' @param rho2 weight placed on the second data set
//' @param lambda the level of sparsity penalization
//' @param eta the level of penalization towards model similarity
//' @param w1 group-specific weights for different penalization across groups in data set 1
//' @param w2 group-specific weights for different penalization across groups in data set 2
//' @param w group-specific weights for different penalization toward similarity for different groups
//' @param AA1 a list of the matrices A1j
//' @param AA1 a list of the matrices A2j
//' @param eigen1 a list of eigen info on groups from data set 1
//' @param eigen2 a list of eigen info on groups from data set 2
//' @param Com the indices of the covariate groups which are common in the two datasets
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta1_init optional starting value for beta1
//' @param beta2_init optional starting value for beta2
//' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
//'
//' @examples
//' data <- get_grouplasso2pop_data(n1 = 400,n2 = 600, response = "binary")
//' 
//' grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = data$Y1,
//'                                                    rX1 = data$X1,
//'                                                    groups1 = data$groups1,
//'                                                    rY2 = data$Y2,
//'                                                    rX2 = data$X2,
//'                                                    groups2 = data$groups2,
//'                                                    rho1 = 2,
//'                                                    rho2 = 1,
//'                                                    lambda = 1,
//'                                                    eta = 1,
//'                                                    w1 = data$w1,
//'                                                    w2 = data$w2,
//'                                                    w = data$w,
//'                                                    rAA1 = data$AA1,
//'                                                    rAA2 = data$AA2,
//'                                                    rCom = data$Com,
//'                                                    tol = 1e-4,
//'                                                    maxiter = 500)
// [[Rcpp::export]]
List grouplasso2pop_logreg_slower(NumericVector rY1,
                                  NumericMatrix rX1,
                                  IntegerVector groups1,
                                  NumericVector rY2,
                                  NumericMatrix rX2,
                                  IntegerVector groups2,
                                  float rho1,
                                  float rho2,
                                  float lambda,
                                  float eta,
                                  NumericVector w1,
                                  NumericVector w2,
                                  NumericVector w,
                                  List rAA1,
                                  List rAA2,
                                  IntegerVector rCom,
                                  float tol,
                                  int maxiter,
                                  NumericVector beta1_init = NumericVector::create(),
                                  NumericVector beta2_init = NumericVector::create()){
  
  // get dimensions
  int n1 = rY1.size(), n2 = rY2.size();
  int q1 = rX1.ncol(), q2 = rX2.ncol();
  
  int i,j,k,l,m;
  
  int g1 = max(groups1);
  int g2 = max(groups2);
  
  LogicalVector got_eigen1(g1);
  LogicalVector got_eigen2(g2);
  
  IntegerVector d1(g1);
  IntegerVector d2(g2);
  IntegerVector is_grp_singleton1(g1);
  IntegerVector is_grp_singleton2(g2);
  IntegerVector grp_begin_ind1(g1);
  IntegerVector grp_begin_ind2(g2);
  
  // make vector indicating what groups are common
  IntegerVector is_com(std::max(g1,g2));
  
  for(j = 0; j < rCom.size(); j++){
    
    is_com[rCom[j]-1] = 1;
    
  }
  
  // identify singleton and non-singleton groups for data set 1 
  j = 1;
  for(i = 0; i < q1; i++){
    
    if(groups1[i] != j){
      
      if( d1[j-1] == 1){
        
        is_grp_singleton1[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d1[j-1]++;
    
  }
  
  grp_begin_ind1[0] = 0;
  for(j = 1; j < g1;j++){
    
    grp_begin_ind1[j] = sum(d1[Range(0,j-1)]);
    
  }
  
  
  // identify singleton and non-singleton groups for data set 2
  j = 1;
  for(i = 0; i < q2; i++){
    
    if(groups2[i] != j){
      
      if( d2[j-1] == 1){
        
        is_grp_singleton2[j-1] = 1;
        
      }
      
      j++;
      
    }
    
    d2[j-1] ++;
    
  }
  
  grp_begin_ind2[0] = 0;
  for(j = 1; j < g2;j++){
    
    grp_begin_ind2[j] = sum(d2[Range(0,j-1)]);
    
  }
  
  // make Armadillo matrices/vectors
  arma::colvec Y1(rY1.begin(),n1,false);
  arma::colvec Y2(rY2.begin(),n2,false);
  
  arma::colvec P1(n1);
  arma::colvec P2(n2);
  
  arma::mat X1(rX1.begin(),n1,q1,false);
  arma::mat X2(rX2.begin(),n2,q2,false);
  
  arma::colvec b1 = arma::zeros(q1);
  arma::colvec b2 = arma::zeros(q2);
  
  // Set up initial values if provided
  if( beta1_init.size() == q1 ){
    
    arma::colvec b1_silly_copy(beta1_init.begin(),q1,false);
    b1 = b1_silly_copy;
    
  }
  
  if( beta2_init.size() == q2 ){
    
    arma::colvec b2_silly_copy(beta2_init.begin(),q2,false);
    b2 = b2_silly_copy;
    
  }
  
  arma::colvec b1_0 = b1;
  arma::colvec b2_0 = b2;
  
  arma::colvec Z1j_tilde(n1);
  arma::colvec h1(n1);
  arma::colvec Z2j_tilde(n2);
  arma::colvec h2(n2);
  
  arma::colvec x1z1j_AAb2(max(d1));
  arma::colvec x2z2j_AAb1(max(d2));
  arma::mat LtL1(max(d1),max(d1));
  arma::mat LtL2(max(d2),max(d2));
  
  arma::field<arma::vec> eigval1(g1);
  arma::field<arma::mat> eigvec1(g1);
  arma::field<arma::mat> cholLtL1(g1);
  
  arma::field<arma::vec> eigval2(g2);
  arma::field<arma::mat> eigvec2(g2);
  arma::field<arma::mat> cholLtL2(g2);
  
  k = max(rCom);
  arma::field<arma::mat> AA1(k);
  arma::field<arma::mat> AA2(k);
  
  for( j = 0; j < k; j++){
    if(is_com[j] == 1){
      
      arma::mat AA1_silly_copy(as<NumericMatrix>(rAA1[j]).begin(),as<NumericMatrix>(rAA1[j]).nrow(),d1[j]);
      AA1(j) = AA1_silly_copy;
      
      arma::mat AA2_silly_copy(as<NumericMatrix>(rAA2[j]).begin(),as<NumericMatrix>(rAA2[j]).nrow(),d2[j]);
      AA2(j) = AA2_silly_copy;
      
    }
    
  }
  
  // define other floats
  float sx1j, sx1z1j_scalar, x1z1j_norm;
  float sx1j_AA, sx1z1j_AAb2_scalar, x1z1j_AAb2_norm;
  
  float sx2j, sx2z2j_scalar, x2z2j_norm;
  float sx2j_AA, sx2z2j_AAb1_scalar, x2z2j_AAb1_norm;
  
  float lambda1 = lambda * n1 / rho1;
  float eta1 = eta * n1 / rho1;
    
  float lambda2 = lambda * n2 / rho1;
  float eta2 = eta * n2 / rho2;
  
  // algorithmic control
  bool conv = false;
  int iter = 0;
  NumericVector obj_val(maxiter);
  
  // begin looping!
  while( (conv == false) & (iter < maxiter)){
    
    b1_0 = b1;
    b2_0 = b2;
    
    // go through groups of data set 1
    for( j = 0; j < g1 ; j++){
      
      // update probs
      P1 = 1/(1 + exp( - X1 * b1  ));
      
      if(any(P1 > 1-1e-10) | any(P1 < 1e-10)){
        
        Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
        iter = maxiter;
        
      }
      
      // first and last columns of X1 belonging to group
      i = grp_begin_ind1[j];
      k = i + d1[j] - 1;
      
      if(is_com[j] == 1){
        
        // first and last columns of X2 belonging to corresponding group
        l = grp_begin_ind2[j];
        m = l + d2[j] - 1;
        
        if(d1[j] == 1){
          
          Z1j_tilde = X1.cols(i,i) * b1(i) / 2 + 2 * ( Y1 - P1);
          sx1j_AA = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i) / 4  + eta1 * w[j] * AA1(j).t() * AA1(j));
          sx1z1j_AAb2_scalar = arma::as_scalar(X1.cols(i,i).t() * Z1j_tilde / 2 + eta1 * w[j]*AA1(j).t() * AA2(j) * b2.rows(l,m));
          
          b1(i) = SoftThresh_scalar(sx1z1j_AAb2_scalar,lambda1*w1[j]/2) / sx1j_AA;  
          
        } else {
          
          Z1j_tilde = X1.cols(i,k) * b1.rows(i,k) / 2 + 2 * ( Y1 - P1);
          
          x1z1j_AAb2.rows(0,d1[j]-1) = X1.cols(i,k).t() * Z1j_tilde / 2 + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m);
          x1z1j_AAb2_norm = sqrt(accu(pow(x1z1j_AAb2.rows(0,d1[j]-1),2)));
          
          if(x1z1j_AAb2_norm < lambda1 * w1[j]/2){
            
            b1.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen1[j] == false){
              
              LtL1.submat(0,0,d1[j]-1,d1[j]-1) = X1.cols(i,k).t() * X1.cols(i,k) / 4 + eta1 * w[j] * AA1(j).t() * AA1(j);
              arma::eig_sym(eigval1(j),eigvec1(j),LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              cholLtL1(j) = chol(LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              got_eigen1[j] = true;
              
            } 
            
            h1 = arma::inv(cholLtL1(j).t()) * x1z1j_AAb2.rows(0,d1[j]-1);
            
            b1.rows(i,k) = FoygelDrton_Armadillo(h1, cholLtL1(j), lambda1 * w1[j]/2, eigval1(j), eigvec1(j).t());
            
          }
          
        }
        
      } else {
        
        if(d1[j] == 1){
          
          
          Z1j_tilde = X1.cols(i,i) * b1(i) / 2 + 2 * ( Y1 - P1);
          sx1j = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i)) / 4;
          sx1z1j_scalar = arma::as_scalar(X1.cols(i,i).t() * Z1j_tilde) / 2;
          b1(i) = SoftThresh_scalar(sx1z1j_scalar,lambda1 * w1[j] / 2) / sx1j;
          
        } else {
          
          Z1j_tilde = X1.cols(i,k) * b1.rows(i,k) / 2 + 2 * ( Y1 - P1);
          x1z1j_norm = sqrt(arma::accu(pow(trans(X1.cols(i,k)) * Z1j_tilde,2))) / 2;
          
          if( x1z1j_norm < lambda1 * w1[j]/2){
            
            b1.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen1[j] == false){
              
              arma::eig_sym(eigval1(j),eigvec1(j),X1.cols(i,k).t() * X1.cols(i,k) / 4);
              got_eigen1[j] = true;
              
            } 
            
            b1.rows(i,k) = FoygelDrton_Armadillo(Z1j_tilde, X1.cols(i,k) / 2, lambda1 * w1[j]/2, eigval1(j), eigvec1(j).t());
            
          }
          
        }
        
      }
      
    }
    
    
    
    // go through groups of data set 2
    for( j = 0; j < g2 ; j++){
      
      // update probs
      P2 = 1/(1 + exp( - X2 * b2  ));
      
      if(any(P2 > 1-1e-10) | any(P2 < 1e-10)){
        
        Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
        iter = maxiter;
        
      }
      
      // first and last columns of X2 belonging to group
      i = grp_begin_ind2[j];
      k = i + d2[j] - 1;
      
      if(is_com[j] == 1){
        
        // first and last columns of X1 belonging to corresponding group
        l = grp_begin_ind1[j];
        m = l + d1[j] - 1;
        
        if(d2[j] == 1){
          
          Z2j_tilde = X2.cols(i,i) * b2(i) / 2 + 2 * ( Y2 - P2);
          sx2j_AA = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i) / 4  + eta2 * w[j] * AA2(j).t() * AA2(j));
          sx2z2j_AAb1_scalar = arma::as_scalar(X2.cols(i,i).t() * Z2j_tilde / 2 + eta2 * w[j]*AA2(j).t() * AA1(j) * b1.rows(l,m));
          
          b2(i) = SoftThresh_scalar(sx2z2j_AAb1_scalar,lambda2*w2[j]/2) / sx2j_AA;  
          
        } else {
          
          Z2j_tilde = X2.cols(i,k) * b2.rows(i,k) / 2 + 2 * ( Y2 - P2);
          
          x2z2j_AAb1.rows(0,d2[j]-1) = X2.cols(i,k).t() * Z2j_tilde / 2 + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m);
          x2z2j_AAb1_norm = sqrt(accu(pow(x2z2j_AAb1.rows(0,d2[j]-1),2)));
          
          if(x2z2j_AAb1_norm < lambda2 * w2[j]/2){
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              LtL2.submat(0,0,d2[j]-1,d2[j]-1) = X2.cols(i,k).t() * X2.cols(i,k) / 4 + eta2 * w[j] * AA2(j).t() * AA2(j);
              arma::eig_sym(eigval2(j),eigvec2(j),LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              cholLtL2(j) = chol(LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              got_eigen2[j] = true;
              
            } 
            
            h2 = arma::inv(cholLtL2(j).t()) * x2z2j_AAb1.rows(0,d2[j]-1);
            
            b2.rows(i,k) = FoygelDrton_Armadillo(h2, cholLtL2(j), lambda2 * w2[j]/2, eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      } else {
        
        if(d2[j] == 1){
          
          Z2j_tilde = X2.cols(i,i) * b2(i) / 2 + 2 * (Y2 - P2);
          sx2j = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i)) / 4;
          sx2z2j_scalar = arma::as_scalar(X2.cols(i,i).t() * Z2j_tilde) / 2;
          b2(i) = SoftThresh_scalar(sx2z2j_scalar,lambda2 * w2[j] / 2) / sx2j;
          
        } else {
          
          Z2j_tilde = X2.cols(i,k) * b2.rows(i,k) / 2 + 2 * (Y2 - P2);
          x2z2j_norm = sqrt(arma::accu(pow(X2.cols(i,k).t() * Z2j_tilde,2))) / 2;
          
          if( x2z2j_norm < lambda2 * w2[j] / 2){
            
            b2.rows(i,k) = arma::zeros(k-i+1);
            
          } else {
            
            if(got_eigen2[j] == false){
              
              arma::eig_sym(eigval2(j),eigvec2(j),X2.cols(i,k).t() * X2.cols(i,k) / 4);
              got_eigen2[j] = true;
              
            } 
            
            b2.rows(i,k) = FoygelDrton_Armadillo(Z2j_tilde, X2.cols(i,k) / 2, lambda2 * w2[j]/2, eigval2(j), eigvec2(j).t());
            
          }
          
        }
        
      }
      
    }
    
    if(any(abs(b1 - b1_0) > tol) | any(abs(b2 - b2_0) > tol) ){
      
      conv = false;
      
    } else {
      
      conv = true;
    }
    
    iter++;
    
  }
  // close while statement
  
  return Rcpp::List::create(Named("beta1.hat") = b1,
                            Named("beta2.hat") = b2,
                            Named("iter") = iter,
                            Named("conv") = conv);
  
}


//' Minimize the objective function of the 2-population group lasso problem with a binary response
//'
//' @param Y1 the binary response vector of data set 1
//' @param X1 matrix containing the design matrices for data set 1
//' @param groups1 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param Y2 the binary response vector of data set 2
//' @param X2 matrix containing the design matrices for data set 2
//' @param groups2 a vector of integers indicating to which group each covariate in data set 1 belongs
//' @param rho1 weight placed on the first data set
//' @param rho2 weight placed on the second data set
//' @param lambda the level of sparsity penalization
//' @param eta the level of penalization towards model similarity
//' @param w1 group-specific weights for different penalization across groups in data set 1
//' @param w2 group-specific weights for different penalization across groups in data set 2
//' @param w group-specific weights for different penalization toward similarity for different groups
//' @param AA1 a list of the matrices A1j
//' @param AA1 a list of the matrices A2j
//' @param eigen1 a list of eigen info on groups from data set 1
//' @param eigen2 a list of eigen info on groups from data set 2
//' @param Com the indices of the covariate groups which are common in the two datasets
//' @param tol a convergence criterion
//' @param max.iter the maximum allowed number of iterations
//' @param return_obj a logical indicating whether the value of the objection function should be recorded after every step of the algorithm
//' @param beta1_init optional starting value for beta1
//' @param beta2_init optional starting value for beta2
//' @return Returns the minimizers of the 2-population group lasso objective function for the two data sets.
//'
//' @examples
//' data <- get_grouplasso2pop_data(n1 = 400,n2 = 600, response = "binary")
//'
//' grouplasso2pop_logreg.out <- grouplasso2pop_logreg(rY1 = data$Y1,
//'                                                    rX1 = data$X1,
//'                                                    groups1 = data$groups1,
//'                                                    rY2 = data$Y2,
//'                                                    rX2 = data$X2,
//'                                                    groups2 = data$groups2,
//'                                                    rho1 = 2,
//'                                                    rho2 = 1,
//'                                                    lambda = 1,
//'                                                    eta = 1,
//'                                                    w1 = data$w1,
//'                                                    w2 = data$w2,
//'                                                    w = data$w,
//'                                                    rAA1 = data$AA1,
//'                                                    rAA2 = data$AA2,
//'                                                    rCom = data$Com,
//'                                                    tol = 1e-4,
//'                                                    maxiter = 500)
// [[Rcpp::export]]
List grouplasso2pop_logreg(NumericVector rY1,
                           NumericMatrix rX1,
                           IntegerVector groups1,
                           NumericVector rY2,
                           NumericMatrix rX2,
                           IntegerVector groups2,
                           float rho1,
                           float rho2,
                           float lambda,
                           float eta,
                           NumericVector w1,
                           NumericVector w2,
                           NumericVector w,
                           List rAA1,
                           List rAA2,
                           IntegerVector rCom,
                           float tol,
                           int maxiter,
                           NumericVector beta1_init = NumericVector::create(),
                           NumericVector beta2_init = NumericVector::create()){

  // get dimensions
  int n1 = rY1.size(), n2 = rY2.size();
  int q1 = rX1.ncol(), q2 = rX2.ncol();

  int i,j,k,l,m;

  int g1 = max(groups1);
  int g2 = max(groups2);

  LogicalVector got_eigen1(g1);
  LogicalVector got_eigen2(g2);

  IntegerVector d1(g1);
  IntegerVector d2(g2);
  IntegerVector is_grp_singleton1(g1);
  IntegerVector is_grp_singleton2(g2);
  IntegerVector grp_begin_ind1(g1);
  IntegerVector grp_begin_ind2(g2);

  // make vector indicating what groups are common
  IntegerVector is_com(std::max(g1,g2));

  for(j = 0; j < rCom.size(); j++){

    is_com[rCom[j]-1] = 1;

  }

  // identify singleton and non-singleton groups for data set 1
  j = 1;
  for(i = 0; i < q1; i++){

    if(groups1[i] != j){

      if( d1[j-1] == 1){

        is_grp_singleton1[j-1] = 1;

      }

      j++;

    }

    d1[j-1]++;

  }

  grp_begin_ind1[0] = 0;
  for(j = 1; j < g1;j++){

    grp_begin_ind1[j] = sum(d1[Range(0,j-1)]);

  }


  // identify singleton and non-singleton groups for data set 2
  j = 1;
  for(i = 0; i < q2; i++){

    if(groups2[i] != j){

      if( d2[j-1] == 1){

        is_grp_singleton2[j-1] = 1;

      }

      j++;

    }

    d2[j-1] ++;

  }

  grp_begin_ind2[0] = 0;
  for(j = 1; j < g2;j++){

    grp_begin_ind2[j] = sum(d2[Range(0,j-1)]);

  }

  // make Armadillo matrices/vectors
  arma::colvec Y1(rY1.begin(),n1,false);
  arma::colvec Y2(rY2.begin(),n2,false);

  arma::colvec P1(n1);
  arma::colvec P2(n2);

  arma::mat X1(rX1.begin(),n1,q1,false);
  arma::mat X2(rX2.begin(),n2,q2,false);

  arma::colvec b1 = arma::zeros(q1);
  arma::colvec b2 = arma::zeros(q2);

  // Set up initial values if provided
  if( beta1_init.size() == q1 ){

    arma::colvec b1_silly_copy(beta1_init.begin(),q1,false);
    b1 = b1_silly_copy;

  }

  if( beta2_init.size() == q2 ){

    arma::colvec b2_silly_copy(beta2_init.begin(),q2,false);
    b2 = b2_silly_copy;

  }

  arma::colvec log_odds1 = X1 * b1;
  arma::colvec log_odds2 = X2 * b2;

  arma::colvec b1_0 = b1;
  arma::colvec b1_00 = b1;
  arma::colvec b2_0 = b2;
  arma::colvec b2_00 = b2;

  arma::colvec Z1j_tilde(n1);
  arma::colvec h1(n1);
  arma::colvec Z2j_tilde(n2);
  arma::colvec h2(n2);

  arma::colvec x1z1j_AAb2(max(d1));
  arma::colvec x2z2j_AAb1(max(d2));
  arma::mat LtL1(max(d1),max(d1));
  arma::mat LtL2(max(d2),max(d2));

  arma::field<arma::vec> eigval1(g1);
  arma::field<arma::mat> eigvec1(g1);
  arma::field<arma::mat> cholLtL1(g1);

  arma::field<arma::vec> eigval2(g2);
  arma::field<arma::mat> eigvec2(g2);
  arma::field<arma::mat> cholLtL2(g2);

  k = max(rCom);
  arma::field<arma::mat> AA1(k);
  arma::field<arma::mat> AA2(k);

  for( j = 0; j < k; j++){
    if(is_com[j] == 1){

      arma::mat AA1_silly_copy(as<NumericMatrix>(rAA1[j]).begin(),as<NumericMatrix>(rAA1[j]).nrow(),d1[j]);
      AA1(j) = AA1_silly_copy;

      arma::mat AA2_silly_copy(as<NumericMatrix>(rAA2[j]).begin(),as<NumericMatrix>(rAA2[j]).nrow(),d2[j]);
      AA2(j) = AA2_silly_copy;

    }

  }

  // define other floats
  float sx1j, sx1z1j_scalar, x1z1j_norm;
  float sx1j_AA, sx1z1j_AAb2_scalar, x1z1j_AAb2_norm;

  float sx2j, sx2z2j_scalar, x2z2j_norm;
  float sx2j_AA, sx2z2j_AAb1_scalar, x2z2j_AAb1_norm;

  float lambda1 = lambda * n1 / rho1;
  float eta1 = eta * n1 / rho1;

  float lambda2 = lambda * n2 / rho1;
  float eta2 = eta * n2 / rho2;

  // algorithmic control
  bool conv = false;
  bool probs0or1 = false;
  int iter = 0;
  NumericVector obj_val(maxiter);

  // begin looping!
  while( (conv == false) & (iter < maxiter)){

    b1_00 = b1;
    b2_00 = b2;

    // go through groups of data set 1
    for( j = 0; j < g1 ; j++){

      // update probs
      P1 = 1/(1 + exp( - log_odds1  ));
      
      if(any(P1 > 1 - 1e-10) | any(P1 < 1e-10)){
        
        iter = maxiter;
        probs0or1 = true;

      }

      // first and last columns of X1 belonging to group
      i = grp_begin_ind1[j];
      k = i + d1[j] - 1;
      b1_0.rows(i,k) = b1.rows(i,k);
      
      if(is_com[j] == 1){

        // first and last columns of X2 belonging to corresponding group
        l = grp_begin_ind2[j];
        m = l + d2[j] - 1;

        if(d1[j] == 1){

          Z1j_tilde = X1.cols(i,i) * b1(i) / 2 + 2 * ( Y1 - P1);
          sx1j_AA = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i) / 4  + eta1 * w[j] * AA1(j).t() * AA1(j));
          sx1z1j_AAb2_scalar = arma::as_scalar(X1.cols(i,i).t() * Z1j_tilde / 2 + eta1 * w[j]*AA1(j).t() * AA2(j) * b2.rows(l,m));

          b1(i) = SoftThresh_scalar(sx1z1j_AAb2_scalar,lambda1*w1[j]/2) / sx1j_AA;

        } else {

          Z1j_tilde = X1.cols(i,k) * b1.rows(i,k) / 2 + 2 * ( Y1 - P1);

          x1z1j_AAb2.rows(0,d1[j]-1) = X1.cols(i,k).t() * Z1j_tilde / 2 + eta1 * w[j] * AA1(j).t() * AA2(j) * b2.rows(l,m);
          x1z1j_AAb2_norm = sqrt(accu(pow(x1z1j_AAb2.rows(0,d1[j]-1),2)));

          if(x1z1j_AAb2_norm < lambda1 * w1[j]/2){

            b1.rows(i,k) = arma::zeros(k-i+1);

          } else {

            if(got_eigen1[j] == false){

              LtL1.submat(0,0,d1[j]-1,d1[j]-1) = X1.cols(i,k).t() * X1.cols(i,k) / 4 + eta1 * w[j] * AA1(j).t() * AA1(j);
              arma::eig_sym(eigval1(j),eigvec1(j),LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              cholLtL1(j) = chol(LtL1.submat(0,0,d1[j]-1,d1[j]-1));
              got_eigen1[j] = true;

            }

            h1 = arma::inv(cholLtL1(j).t()) * x1z1j_AAb2.rows(0,d1[j]-1);

            b1.rows(i,k) = FoygelDrton_Armadillo(h1, cholLtL1(j), lambda1 * w1[j]/2, eigval1(j), eigvec1(j).t());

          }

        }

      } else {

        if(d1[j] == 1){


          Z1j_tilde = X1.cols(i,i) * b1(i) / 2 + 2 * ( Y1 - P1);
          sx1j = arma::as_scalar(X1.cols(i,i).t() * X1.cols(i,i)) / 4;
          sx1z1j_scalar = arma::as_scalar(X1.cols(i,i).t() * Z1j_tilde) / 2;
          b1(i) = SoftThresh_scalar(sx1z1j_scalar,lambda1 * w1[j] / 2) / sx1j;

        } else {

          Z1j_tilde = X1.cols(i,k) * b1.rows(i,k) / 2 + 2 * ( Y1 - P1);
          x1z1j_norm = sqrt(arma::accu(pow(trans(X1.cols(i,k)) * Z1j_tilde,2))) / 2;

          if( x1z1j_norm < lambda1 * w1[j]/2){

            b1.rows(i,k) = arma::zeros(k-i+1);

          } else {

            if(got_eigen1[j] == false){

              arma::eig_sym(eigval1(j),eigvec1(j),X1.cols(i,k).t() * X1.cols(i,k) / 4);
              got_eigen1[j] = true;

            }

            b1.rows(i,k) = FoygelDrton_Armadillo(Z1j_tilde, X1.cols(i,k) / 2, lambda1 * w1[j]/2, eigval1(j), eigvec1(j).t());

          }

        }

      }

      log_odds1 = log_odds1 - X1.cols(i,k) * (b1_0.rows(i,k) - b1.rows(i,k));
      
    }


    // go through groups of data set 2
    for( j = 0; j < g2 ; j++){

      // update probs
      P2 = 1/(1 + exp( - log_odds2  ));

      if(any(P2 > 1 - 1e-10) | any(P2 < 1e-10)){
        
        iter = maxiter;
        probs0or1 = true;
        
      }
      
      // first and last columns of X2 belonging to group
      i = grp_begin_ind2[j];
      k = i + d2[j] - 1;
      b2_0.rows(i,k) = b2.rows(i,k);

      if(is_com[j] == 1){

        // first and last columns of X1 belonging to corresponding group
        l = grp_begin_ind1[j];
        m = l + d1[j] - 1;

        if(d2[j] == 1){

          Z2j_tilde = X2.cols(i,i) * b2(i) / 2 + 2 * ( Y2 - P2);
          sx2j_AA = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i) / 4  + eta2 * w[j] * AA2(j).t() * AA2(j));
          sx2z2j_AAb1_scalar = arma::as_scalar(X2.cols(i,i).t() * Z2j_tilde / 2 + eta2 * w[j]*AA2(j).t() * AA1(j) * b1.rows(l,m));

          b2(i) = SoftThresh_scalar(sx2z2j_AAb1_scalar,lambda2*w2[j]/2) / sx2j_AA;

        } else {

          Z2j_tilde = X2.cols(i,k) * b2.rows(i,k) / 2 + 2 * ( Y2 - P2);

          x2z2j_AAb1.rows(0,d2[j]-1) = X2.cols(i,k).t() * Z2j_tilde / 2 + eta2 * w[j] * AA2(j).t() * AA1(j) * b1.rows(l,m);
          x2z2j_AAb1_norm = sqrt(accu(pow(x2z2j_AAb1.rows(0,d2[j]-1),2)));

          if(x2z2j_AAb1_norm < lambda2 * w2[j]/2){

            b2.rows(i,k) = arma::zeros(k-i+1);

          } else {

            if(got_eigen2[j] == false){

              LtL2.submat(0,0,d2[j]-1,d2[j]-1) = X2.cols(i,k).t() * X2.cols(i,k) / 4 + eta2 * w[j] * AA2(j).t() * AA2(j);
              arma::eig_sym(eigval2(j),eigvec2(j),LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              cholLtL2(j) = chol(LtL2.submat(0,0,d2[j]-1,d2[j]-1));
              got_eigen2[j] = true;

            }

            h2 = arma::inv(cholLtL2(j).t()) * x2z2j_AAb1.rows(0,d2[j]-1);

            b2.rows(i,k) = FoygelDrton_Armadillo(h2, cholLtL2(j), lambda2 * w2[j]/2, eigval2(j), eigvec2(j).t());

          }

        }

      } else {

        if(d2[j] == 1){

          Z2j_tilde = X2.cols(i,i) * b2(i) / 2 + 2 * (Y2 - P2);
          sx2j = arma::as_scalar(X2.cols(i,i).t() * X2.cols(i,i)) / 4;
          sx2z2j_scalar = arma::as_scalar(X2.cols(i,i).t() * Z2j_tilde) / 2;
          b2(i) = SoftThresh_scalar(sx2z2j_scalar,lambda2 * w2[j] / 2) / sx2j;

        } else {

          Z2j_tilde = X2.cols(i,k) * b2.rows(i,k) / 2 + 2 * (Y2 - P2);
          x2z2j_norm = sqrt(arma::accu(pow(X2.cols(i,k).t() * Z2j_tilde,2))) / 2;

          if( x2z2j_norm < lambda2 * w2[j] / 2){

            b2.rows(i,k) = arma::zeros(k-i+1);

          } else {

            if(got_eigen2[j] == false){

              arma::eig_sym(eigval2(j),eigvec2(j),X2.cols(i,k).t() * X2.cols(i,k) / 4);
              got_eigen2[j] = true;

            }

            b2.rows(i,k) = FoygelDrton_Armadillo(Z2j_tilde, X2.cols(i,k) / 2, lambda2 * w2[j]/2, eigval2(j), eigvec2(j).t());

          }

        }

      }
      
      log_odds2 = log_odds2 - X2.cols(i,k) * (b2_0.rows(i,k) - b2.rows(i,k));
      
    }

    if(any(abs(b1 - b1_00) > tol) | any(abs(b2 - b2_00) > tol) ){

      conv = false;

    } else {

      conv = true;
      
    }

    iter++;

  }
  
  // close while statement
  
  if(probs0or1){
    
    Rcpp::Rcout << "warning: fitted probabilities tending to 0 or 1" << std::endl;
    
  }
  

  return Rcpp::List::create(Named("beta1.hat") = b1,
                            Named("beta2.hat") = b2,
                            Named("iter") = iter,
                            Named("conv") = conv,
                            Named("probs0or1") = probs0or1);

}

//' Generate all possible sequences of 0s and 1s of a given length
//' @param a the length of the sequences.
//' @return a matrix containing in its rows the sequences of 0s and 1s.
// [[Rcpp::export]]
arma::mat all_binary_sequences(int a){
  
  int n_rows = (int)(pow(2,a) + .5);
  int i,n,m;
  arma::mat M = arma::zeros(n_rows,a);
  
  // count from 0 to 2^a - 1 in binary
  for( n = 0; n < n_rows ; n++ ){
    
    m = n;
    i = 0;
    while( m > 0){
      
      M(n,i) = m % 2;    
      m = m / 2;  
      i++;
      
    }  
    
  }
  
  return M;
  
}

//' Computes conditional expectations of individual disease statuses for individual, master pool, or Dorfman testing
//'   
//' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
//' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
//' @param X Design matrix with first column a column of 1s.
//' @param b Parameter values at which to compute the conditional expectations.
//' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
//' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
//' @return The vector of conditional expectations.
//' 
//' This function computes the conditional expectations of each individual disease status, conditional on the observed assay data and the diseasestatuses of all other individuals.
//' 
//' @examples
//' grouplassogt2pop_data <- get_grouplassogt2pop_data( n1 = 400, n2 = 600)
//'   
//' EY <- EYexact(Z = grouplassogt2pop_data$Z1,
//'               Y = grouplassogt2pop_data$Y1,
//'               X = grouplassogt2pop_data$X1,
//'               b = rep(1,ncol(grouplassogt2pop_data$X1)),
//'               Se = grouplassogt2pop_data$Se1,
//'               Sp = grouplassogt2pop_data$Sp1)
// [[Rcpp::export]]
arma::colvec EYexact(IntegerMatrix Z, IntegerMatrix Y, NumericMatrix X, NumericVector b, NumericVector Se, NumericVector Sp){
  
  int n = X.nrow(), p = X.ncol(), max_c = Z.ncol() - 3;
  int assay_in_which_involved;
  int group_assay;
  int i, j, k, l, cj, zj;
  int which_SeSp;
  int n_left,n_group;
  
  double Sej, Spj, prd, max_Y;
  double prob_whole_group_negative;
  
  // arma versions of inputs
  arma::imat aY(Y.begin(),n,4,false);
  arma::imat aZ(Z.begin(),Z.nrow(),Z.ncol(),false);
  
  // initialize some armadillo vectors
  arma::ivec group_long(max_c);
  arma::ivec U(max_c);
  arma::uvec individual_assays(max_c);
  arma::uvec group(max_c);
  arma::uvec to_drop(1);
  arma::colvec Sej_i(max_c);
  arma::colvec Spj_i(max_c);
  arma::colvec p_group(max_c);
  arma::colvec Aj(max_c);
  arma::colvec Bj(max_c);
  arma::colvec EY(n);
  arma::colvec aa(max_c);
  arma::colvec bb(max_c);
  arma::colvec cc(max_c);
  arma::colvec aa_bb_cc(max_c);
  arma::colvec Y_no_i(max_c-1);
  arma::umat cj_by_cj(max_c-1,max_c);
  arma::uvec k_uvec(1);
  
  // get the fitted probs
  arma::mat aX(X.begin(),n,p,false);
  arma::colvec ab(b.begin(),p,false);
  arma::colvec p_all = 1 / ( 1 + exp( - aX * ab )) ;
  
  // make lists of individuals tested in a single assay and in two assays
  arma::uvec involved_in_only_one_assay = arma::find( aY.col(1) == 1 );
  arma::uvec involved_in_two_assays = arma::find( aY.col(1) == 2 );
  
  // generate matrices with rows giving all sequences of 0s and 1s of lengths cj
  arma::field<arma::mat> YY(max_c - 1);
  for(i = 0; i < max_c - 1; i++){
    
    YY(i) = all_binary_sequences(i + 2);
    
  }
  
  //----------------------------------------------------------
  // handle individual, masterpool, and negative dorfman pools
  //----------------------------------------------------------
  
  n_left = involved_in_only_one_assay.size();
  while(n_left > 0){
    
    // in which assay was he/she involved?
    assay_in_which_involved = aY(involved_in_only_one_assay[0],2) - 1;
    
    // all involved in this assay
    group_long = aZ(assay_in_which_involved,arma::span(3,3 + max_c - 1)).t() - 1;
    group.set_size(max_c);
    cj = 0;
    for(i = 0; i < max_c ; i++){
      
      if(group_long(i) == -100) break;
      
      group(i) = group_long(i);
      cj++;
      
    }
    
    // set size of some arma vectors
    group.set_size(cj);
    p_group.set_size(cj);
    Aj.set_size(cj);
    Bj.set_size(cj);
    
    // get probs for the group
    p_group = p_all(group);
    
    // get Se and Sp for this assay
    which_SeSp = aZ(assay_in_which_involved,2) - 1;
    Sej = Se[which_SeSp];
    Spj = Sp[which_SeSp];
    
    // prepare to compute the group probabilities
    prob_whole_group_negative = prod(1 - p_group);
    zj = aZ(assay_in_which_involved,0);
    Aj = arma::ones(cj,1) * (Sej*zj + (1-Sej)*(1-zj));
    Bj = ((1-Spj)*zj + Spj*(1-zj) - Aj ) * prob_whole_group_negative / (1 - p_group) + Aj;
    
    // compute group probabilities
    EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
    
    // remove the "processed" individuals from the waiting list
    n_group = group.size();
    
    for(i = 0; i < n_group; i++){
      
      to_drop = find(involved_in_only_one_assay == group(i));
      involved_in_only_one_assay.shed_row(to_drop[0]);
      
    }
    
    n_left = involved_in_only_one_assay.size();
    
  }
  
  //--------------------------------------------
  // handle positive pools under dorfman testing
  //--------------------------------------------
  
  n_left = involved_in_two_assays.size();
  while(n_left > 0){
    
    // in which assays was he/she involved? ASSUMING THE GROUP ASSAY COMES FIRST!!! 
    group_assay = aY(involved_in_two_assays[0],2) - 1;
    
    // all involved in this group : ASSUMING THE GROUP ASSAY COMES FIRST!!! 
    // Which it should, unless the data are encoded in the Z matrix in a different way.
    group_long = aZ(group_assay,arma::span(3,3 + max_c - 1)).t() - 1;
    group.set_size(max_c);
    cj = 0;
    for(i = 0; i < max_c ; i++){
      
      if(group_long(i) == -100) break;
      
      group(i) = group_long(i);
      cj++;
      
    }
    
    // reset size of some arma vectors/matrices
    group.set_size(cj);
    p_group.set_size(cj);
    Aj.set_size(cj);
    Bj.set_size(cj);
    individual_assays.set_size(cj);
    Sej_i.set_size(cj);
    Spj_i.set_size(cj);
    U.set_size(cj);
    aa.set_size(cj);
    bb.set_size(cj);
    cc.set_size(cj);
    aa_bb_cc.set_size(cj);
    cj_by_cj.set_size(cj-1,cj);
    U.set_size(cj);
    
    // get probs for the group
    p_group = p_all(group);
    
    // get Se and Sp for the group assay:
    which_SeSp = aZ(group_assay,2) - 1;
    Sej = Se[which_SeSp];
    Spj = Sp[which_SeSp];
    
    // get individual assays and Se and Sp for individual assays
    for( i = 0; i < cj ; i ++){
      
      individual_assays(i) = aY(group(i),3) - 1;
      which_SeSp = aZ(individual_assays(i),2) - 1;
      Sej_i(i) = Se[which_SeSp];
      Spj_i(i) = Sp[which_SeSp];
      U(i) = aZ(individual_assays(i),0);
      
    }
    
    Aj = arma::zeros(cj);
    Bj = arma::zeros(cj);
    
    for(k = 0; k < YY(cj-2).n_rows ; k++){
      
      k_uvec(0) = k;
      
      int yi;
      
      for(i = 0; i < cj ; i++){
        
        yi = YY(cj-2)(k,i);
        
        if( yi == 1){
          
          aa(i) = (Sej_i(i) * U(i) + ( 1 - Sej_i(i)) * (1-U(i)));
          bb(i) = 1;
          cc(i) = p_group(i);
          
        } else {
          
          aa(i) = 1;
          bb(i) = ((1-Spj_i(i))*U(i) + Spj_i(i)*(1-U(i)));
          cc(i) = 1 - p_group(i);
        }
        
        l = 0;
        j = 0;
        while(l < cj - 1){
          
          cj_by_cj(l,i) = j;
          
          if(j != i){
            l++;
          }
          
          j++;
          
        }
        
      }
      
      aa_bb_cc = aa % bb % cc;
      
      for(i = 0; i < cj; i++){
        
        prd = prod(aa_bb_cc(cj_by_cj.col(i)));
        Y_no_i = YY(cj-2)(k_uvec,cj_by_cj.col(i)).t();
        max_Y = max(Y_no_i);
        
        Aj(i) = Aj(i) + prd;
        Bj(i) = Bj(i) + (Sej*max_Y + (1-Spj)*(1 - max_Y)) * prd;
        
      }
      
    }
    
    Aj = Sej * (Sej_i % U + ( 1 - Sej_i) % (1 - U) ) % Aj;
    Bj = ( (1-Spj_i) % U + Spj_i % (1-U) ) % Bj;
    
    EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
    
    // remove the "processed" individuals from the waiting list
    n_group = group.size();
    
    for(i = 0; i < n_group; i++){
      
      to_drop = find(involved_in_two_assays == group(i));
      involved_in_two_assays.shed_row(to_drop[0]);
      
    }
    
    n_left = involved_in_two_assays.size();
    
  }
  
  return EY;
  
}
//' Computes approximate conditional expectations of individual disease statuses for individual, master pool, or Dorfman testing
//'   
//' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
//' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
//' @param X Design matrix with first column a column of 1s.
//' @param b Parameter values at which to compute the conditional expectations.
//' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
//' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
//' @return The vector of conditional expectations.
//' 
//' This function computes approximate conditional expectations via Gibbs sampling of each individual disease status, conditional on the observed assay data and the diseasestatuses of all other individuals.
//' 
//' @examples
//' grouplassogt2pop_data <- get_grouplassogt2pop_data( n1 = 400, n2 = 600)
//'   
//' EY <- EYapprox(Z = grouplassogt2pop_data$Z1,
//'                Y = grouplassogt2pop_data$Y1,
//'                X = grouplassogt2pop_data$X1,
//'                b = rep(1,ncol(grouplassogt2pop_data$X1)),
//'                Se = grouplassogt2pop_data$Se1,
//'                Sp = grouplassogt2pop_data$Sp1)
// [[Rcpp::export]]
Rcpp::NumericVector EYgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z,
                            NumericVector se, NumericVector sp, 
                            int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj, u;
  NumericVector WW(N);
  
  for(g=0;g<GI;g++){
    for(i=0;i<N;i++){
      pi1=p(i);
      pi2=1-p(i);
      np=Y(i,1);
      for(l=0;l<np;l++){
        j=Y(i,(l+2));
        Zj=Z(j-1,0);
        cj=Z(j-1,1);
        tid=Z(j-1,2);
        sej=se(tid-1);
        spj=sp(tid-1);
        ybar=0;
        Y(i,0)=0;
        for(t=0;t<cj;t++){
          id=Z(j-1,(t+3));
          ybar=ybar+Y(id-1,0);
        }
        pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
        if(ybar > 0){
          pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
        }else{
          pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
        }
      }
      pistar=(pi1/(pi1+pi2));
      u = R::runif(0,1);
      //    u=rand() / double(RAND_MAX);
      if(u<pistar){
        Y(i,0)=1;
      }else{Y(i,0)=0;}
      WW(i)=WW(i)+Y(i,0);
    }}  
  
  return WW;
}





