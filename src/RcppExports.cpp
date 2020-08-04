// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SoftThresh_scalar
float SoftThresh_scalar(float z, float a);
RcppExport SEXP _semipaddgt2pop_SoftThresh_scalar(SEXP zSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type z(zSEXP);
    Rcpp::traits::input_parameter< float >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(SoftThresh_scalar(z, a));
    return rcpp_result_gen;
END_RCPP
}
// FoygelDrton_Armadillo
arma::colvec FoygelDrton_Armadillo(arma::colvec h, arma::mat L, double lambda, arma::colvec evals, arma::mat evecs);
RcppExport SEXP _semipaddgt2pop_FoygelDrton_Armadillo(SEXP hSEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP evalsSEXP, SEXP evecsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type evals(evalsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type evecs(evecsSEXP);
    rcpp_result_gen = Rcpp::wrap(FoygelDrton_Armadillo(h, L, lambda, evals, evecs));
    return rcpp_result_gen;
END_RCPP
}
// grouplasso_linreg
List grouplasso_linreg(NumericVector rY, NumericMatrix rX, IntegerVector groups, float lambda, NumericVector w, float tol, int maxiter, NumericVector beta_init);
RcppExport SEXP _semipaddgt2pop_grouplasso_linreg(SEXP rYSEXP, SEXP rXSEXP, SEXP groupsSEXP, SEXP lambdaSEXP, SEXP wSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP beta_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rY(rYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX(rXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_init(beta_initSEXP);
    rcpp_result_gen = Rcpp::wrap(grouplasso_linreg(rY, rX, groups, lambda, w, tol, maxiter, beta_init));
    return rcpp_result_gen;
END_RCPP
}
// grouplasso2pop_linreg
List grouplasso2pop_linreg(NumericVector rY1, NumericMatrix rX1, IntegerVector groups1, NumericVector rY2, NumericMatrix rX2, IntegerVector groups2, float rho1, float rho2, float lambda, float eta, NumericVector w1, NumericVector w2, NumericVector w, List rAA1, List rAA2, IntegerVector rCom, float tol, int maxiter, NumericVector beta1_init, NumericVector beta2_init);
RcppExport SEXP _semipaddgt2pop_grouplasso2pop_linreg(SEXP rY1SEXP, SEXP rX1SEXP, SEXP groups1SEXP, SEXP rY2SEXP, SEXP rX2SEXP, SEXP groups2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP lambdaSEXP, SEXP etaSEXP, SEXP w1SEXP, SEXP w2SEXP, SEXP wSEXP, SEXP rAA1SEXP, SEXP rAA2SEXP, SEXP rComSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP beta1_initSEXP, SEXP beta2_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rY1(rY1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX1(rX1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups1(groups1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rY2(rY2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX2(rX2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups2(groups2SEXP);
    Rcpp::traits::input_parameter< float >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< float >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< float >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type rAA1(rAA1SEXP);
    Rcpp::traits::input_parameter< List >::type rAA2(rAA2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rCom(rComSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta1_init(beta1_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta2_init(beta2_initSEXP);
    rcpp_result_gen = Rcpp::wrap(grouplasso2pop_linreg(rY1, rX1, groups1, rY2, rX2, groups2, rho1, rho2, lambda, eta, w1, w2, w, rAA1, rAA2, rCom, tol, maxiter, beta1_init, beta2_init));
    return rcpp_result_gen;
END_RCPP
}
// grouplasso_logreg
List grouplasso_logreg(NumericVector rY, NumericMatrix rX, IntegerVector groups, float lambda, NumericVector w, float tol, int maxiter, NumericVector beta_init);
RcppExport SEXP _semipaddgt2pop_grouplasso_logreg(SEXP rYSEXP, SEXP rXSEXP, SEXP groupsSEXP, SEXP lambdaSEXP, SEXP wSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP beta_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rY(rYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX(rXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_init(beta_initSEXP);
    rcpp_result_gen = Rcpp::wrap(grouplasso_logreg(rY, rX, groups, lambda, w, tol, maxiter, beta_init));
    return rcpp_result_gen;
END_RCPP
}
// grouplasso2pop_logreg
List grouplasso2pop_logreg(NumericVector rY1, NumericMatrix rX1, IntegerVector groups1, NumericVector rY2, NumericMatrix rX2, IntegerVector groups2, float rho1, float rho2, float lambda, float eta, NumericVector w1, NumericVector w2, NumericVector w, List rAA1, List rAA2, IntegerVector rCom, float tol, int maxiter, NumericVector beta1_init, NumericVector beta2_init);
RcppExport SEXP _semipaddgt2pop_grouplasso2pop_logreg(SEXP rY1SEXP, SEXP rX1SEXP, SEXP groups1SEXP, SEXP rY2SEXP, SEXP rX2SEXP, SEXP groups2SEXP, SEXP rho1SEXP, SEXP rho2SEXP, SEXP lambdaSEXP, SEXP etaSEXP, SEXP w1SEXP, SEXP w2SEXP, SEXP wSEXP, SEXP rAA1SEXP, SEXP rAA2SEXP, SEXP rComSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP beta1_initSEXP, SEXP beta2_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rY1(rY1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX1(rX1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups1(groups1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rY2(rY2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rX2(rX2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groups2(groups2SEXP);
    Rcpp::traits::input_parameter< float >::type rho1(rho1SEXP);
    Rcpp::traits::input_parameter< float >::type rho2(rho2SEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< float >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type rAA1(rAA1SEXP);
    Rcpp::traits::input_parameter< List >::type rAA2(rAA2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rCom(rComSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta1_init(beta1_initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta2_init(beta2_initSEXP);
    rcpp_result_gen = Rcpp::wrap(grouplasso2pop_logreg(rY1, rX1, groups1, rY2, rX2, groups2, rho1, rho2, lambda, eta, w1, w2, w, rAA1, rAA2, rCom, tol, maxiter, beta1_init, beta2_init));
    return rcpp_result_gen;
END_RCPP
}
// all_binary_sequences
arma::mat all_binary_sequences(int a);
RcppExport SEXP _semipaddgt2pop_all_binary_sequences(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(all_binary_sequences(a));
    return rcpp_result_gen;
END_RCPP
}
// EYexact
arma::colvec EYexact(IntegerMatrix Z, IntegerMatrix Y, NumericMatrix X, NumericVector b, NumericVector Se, NumericVector Sp);
RcppExport SEXP _semipaddgt2pop_EYexact(SEXP ZSEXP, SEXP YSEXP, SEXP XSEXP, SEXP bSEXP, SEXP SeSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Se(SeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(EYexact(Z, Y, X, b, Se, Sp));
    return rcpp_result_gen;
END_RCPP
}
// EYgibbs
Rcpp::NumericVector EYgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericVector se, NumericVector sp, int na, int GI);
RcppExport SEXP _semipaddgt2pop_EYgibbs(SEXP NSEXP, SEXP pSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP seSEXP, SEXP spSEXP, SEXP naSEXP, SEXP GISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se(seSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< int >::type na(naSEXP);
    Rcpp::traits::input_parameter< int >::type GI(GISEXP);
    rcpp_result_gen = Rcpp::wrap(EYgibbs(N, p, Y, Z, se, sp, na, GI));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_semipaddgt2pop_SoftThresh_scalar", (DL_FUNC) &_semipaddgt2pop_SoftThresh_scalar, 2},
    {"_semipaddgt2pop_FoygelDrton_Armadillo", (DL_FUNC) &_semipaddgt2pop_FoygelDrton_Armadillo, 5},
    {"_semipaddgt2pop_grouplasso_linreg", (DL_FUNC) &_semipaddgt2pop_grouplasso_linreg, 8},
    {"_semipaddgt2pop_grouplasso2pop_linreg", (DL_FUNC) &_semipaddgt2pop_grouplasso2pop_linreg, 20},
    {"_semipaddgt2pop_grouplasso_logreg", (DL_FUNC) &_semipaddgt2pop_grouplasso_logreg, 8},
    {"_semipaddgt2pop_grouplasso2pop_logreg", (DL_FUNC) &_semipaddgt2pop_grouplasso2pop_logreg, 20},
    {"_semipaddgt2pop_all_binary_sequences", (DL_FUNC) &_semipaddgt2pop_all_binary_sequences, 1},
    {"_semipaddgt2pop_EYexact", (DL_FUNC) &_semipaddgt2pop_EYexact, 6},
    {"_semipaddgt2pop_EYgibbs", (DL_FUNC) &_semipaddgt2pop_EYgibbs, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_semipaddgt2pop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
