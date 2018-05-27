#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DPMMmlogit_chol_rcpp(SEXP);
extern SEXP _DPMMmlogit_invHess(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DPMMmlogit_log_dmvnorm_arma_mult(SEXP, SEXP, SEXP);
extern SEXP _DPMMmlogit_log_dmvnorm_arma_sing(SEXP, SEXP, SEXP);
extern SEXP _DPMMmlogit_logdiwish_rcpp(SEXP, SEXP, SEXP);
extern SEXP _DPMMmlogit_prob_mlogit(SEXP, SEXP);
extern SEXP _DPMMmlogit_riwish_rcpp(SEXP, SEXP);
extern SEXP _DPMMmlogit_riwish_rcpp_chol(SEXP, SEXP);
extern SEXP _DPMMmlogit_rmvnormrcpp(SEXP, SEXP, SEXP);
extern SEXP _DPMMmlogit_rmvnormrcpp_chol(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_DPMMmlogit_chol_rcpp",             (DL_FUNC) &_DPMMmlogit_chol_rcpp,             1},
    {"_DPMMmlogit_invHess",               (DL_FUNC) &_DPMMmlogit_invHess,               5},
    {"_DPMMmlogit_log_dmvnorm_arma_mult", (DL_FUNC) &_DPMMmlogit_log_dmvnorm_arma_mult, 3},
    {"_DPMMmlogit_log_dmvnorm_arma_sing", (DL_FUNC) &_DPMMmlogit_log_dmvnorm_arma_sing, 3},
    {"_DPMMmlogit_logdiwish_rcpp",        (DL_FUNC) &_DPMMmlogit_logdiwish_rcpp,        3},
    {"_DPMMmlogit_prob_mlogit",           (DL_FUNC) &_DPMMmlogit_prob_mlogit,           2},
    {"_DPMMmlogit_riwish_rcpp",           (DL_FUNC) &_DPMMmlogit_riwish_rcpp,           2},
    {"_DPMMmlogit_riwish_rcpp_chol",      (DL_FUNC) &_DPMMmlogit_riwish_rcpp_chol,      2},
    {"_DPMMmlogit_rmvnormrcpp",           (DL_FUNC) &_DPMMmlogit_rmvnormrcpp,           3},
    {"_DPMMmlogit_rmvnormrcpp_chol",      (DL_FUNC) &_DPMMmlogit_rmvnormrcpp_chol,      3},
    {NULL, NULL, 0}
};

void R_init_DPMMmlogit(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
