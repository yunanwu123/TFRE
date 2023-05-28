// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// QICD
Eigen::MatrixXd QICD(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd lambda_list, Eigen::VectorXd initial, const double thresh, const int maxin, const int maxout);
RcppExport SEXP _TFRE_QICD(SEXP XSEXP, SEXP ySEXP, SEXP lambda_listSEXP, SEXP initialSEXP, SEXP threshSEXP, SEXP maxinSEXP, SEXP maxoutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type lambda_list(lambda_listSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< const double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< const int >::type maxin(maxinSEXP);
    Rcpp::traits::input_parameter< const int >::type maxout(maxoutSEXP);
    rcpp_result_gen = Rcpp::wrap(QICD(X, y, lambda_list, initial, thresh, maxin, maxout));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TFRE_QICD", (DL_FUNC) &_TFRE_QICD, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_TFRE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
