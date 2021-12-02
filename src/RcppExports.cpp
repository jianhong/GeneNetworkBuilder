// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// filterNodes
RcppExport SEXP filterNodes(SEXP xx_from, SEXP xx_to, SEXP xx_miRNA, SEXP xx_logFC, SEXP xx_pval, SEXP xx_dir, SEXP rows, SEXP rootgene, SEXP rootlogFC, SEXP tol, SEXP minify, SEXP miRNAtol, SEXP lFC, SEXP pVAL);
RcppExport SEXP _GeneNetworkBuilder_filterNodes(SEXP xx_fromSEXP, SEXP xx_toSEXP, SEXP xx_miRNASEXP, SEXP xx_logFCSEXP, SEXP xx_pvalSEXP, SEXP xx_dirSEXP, SEXP rowsSEXP, SEXP rootgeneSEXP, SEXP rootlogFCSEXP, SEXP tolSEXP, SEXP minifySEXP, SEXP miRNAtolSEXP, SEXP lFCSEXP, SEXP pVALSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xx_from(xx_fromSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xx_to(xx_toSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xx_miRNA(xx_miRNASEXP);
    Rcpp::traits::input_parameter< SEXP >::type xx_logFC(xx_logFCSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xx_pval(xx_pvalSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xx_dir(xx_dirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rootgene(rootgeneSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rootlogFC(rootlogFCSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< SEXP >::type minify(minifySEXP);
    Rcpp::traits::input_parameter< SEXP >::type miRNAtol(miRNAtolSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lFC(lFCSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pVAL(pVALSEXP);
    rcpp_result_gen = Rcpp::wrap(filterNodes(xx_from, xx_to, xx_miRNA, xx_logFC, xx_pval, xx_dir, rows, rootgene, rootlogFC, tol, minify, miRNAtol, lFC, pVAL));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP filterNodes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GeneNetworkBuilder_filterNodes", (DL_FUNC) &_GeneNetworkBuilder_filterNodes, 14},
    {"filterNodes", (DL_FUNC) &filterNodes, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_GeneNetworkBuilder(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
