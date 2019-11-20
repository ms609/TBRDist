// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// uspr_dist
IntegerVector uspr_dist(StringVector tree1, StringVector tree2, LogicalVector keepLabels, LogicalVector useTbrApproxEstimate, LogicalVector useTbrEstimate, LogicalVector useReplugEstimate);
RcppExport SEXP _TBRDist_uspr_dist(SEXP tree1SEXP, SEXP tree2SEXP, SEXP keepLabelsSEXP, SEXP useTbrApproxEstimateSEXP, SEXP useTbrEstimateSEXP, SEXP useReplugEstimateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type tree1(tree1SEXP);
    Rcpp::traits::input_parameter< StringVector >::type tree2(tree2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type keepLabels(keepLabelsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type useTbrApproxEstimate(useTbrApproxEstimateSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type useTbrEstimate(useTbrEstimateSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type useReplugEstimate(useReplugEstimateSEXP);
    rcpp_result_gen = Rcpp::wrap(uspr_dist(tree1, tree2, keepLabels, useTbrApproxEstimate, useTbrEstimate, useReplugEstimate));
    return rcpp_result_gen;
END_RCPP
}
// tbr_dist
List tbr_dist(StringVector tree1, StringVector tree2, LogicalVector printMafs, LogicalVector countMafs, LogicalVector keepLabels, LogicalVector optimize, LogicalVector protectB, LogicalVector exact, LogicalVector approximate);
RcppExport SEXP _TBRDist_tbr_dist(SEXP tree1SEXP, SEXP tree2SEXP, SEXP printMafsSEXP, SEXP countMafsSEXP, SEXP keepLabelsSEXP, SEXP optimizeSEXP, SEXP protectBSEXP, SEXP exactSEXP, SEXP approximateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type tree1(tree1SEXP);
    Rcpp::traits::input_parameter< StringVector >::type tree2(tree2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type printMafs(printMafsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type countMafs(countMafsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type keepLabels(keepLabelsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type optimize(optimizeSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type protectB(protectBSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type approximate(approximateSEXP);
    rcpp_result_gen = Rcpp::wrap(tbr_dist(tree1, tree2, printMafs, countMafs, keepLabels, optimize, protectB, exact, approximate));
    return rcpp_result_gen;
END_RCPP
}
// replug_dist
List replug_dist(StringVector tree1, StringVector tree2, LogicalVector keepLabels);
RcppExport SEXP _TBRDist_replug_dist(SEXP tree1SEXP, SEXP tree2SEXP, SEXP keepLabelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type tree1(tree1SEXP);
    Rcpp::traits::input_parameter< StringVector >::type tree2(tree2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type keepLabels(keepLabelsSEXP);
    rcpp_result_gen = Rcpp::wrap(replug_dist(tree1, tree2, keepLabels));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TBRDist_uspr_dist", (DL_FUNC) &_TBRDist_uspr_dist, 6},
    {"_TBRDist_tbr_dist", (DL_FUNC) &_TBRDist_tbr_dist, 9},
    {"_TBRDist_replug_dist", (DL_FUNC) &_TBRDist_replug_dist, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_TBRDist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}