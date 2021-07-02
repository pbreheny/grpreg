#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
SEXP getListElement(SEXP list, const char *str);
double pchisq(double x, double df, int lower_tail, int give_log);

SEXP mfdr_gaussian(SEXP fit) {
  
  // Declarations
  int n = INTEGER(getListElement(fit, "n"))[0];
  int L = ncols(getListElement(fit, "beta"));
  int ng = size(getListElement(fit, "group.multiplier"));  
  int *gm = nrows(getListElement(fit, "group.multiplier"));  
  double *b = REAL(getListElement(fit, "beta"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double *df = REAL(getListElement(fit, "df"));
  double *RSS = REAL(getListElement(fit, "loss"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double *m = REAL(getListElement(fit, "penalty.factor"));
  double tauSq;
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;
  
  // Calculation
  for (int l=0; l<L; l++) {
    for (int j=1; j < ng; j++) {
      tauSq = RSS[l]/(n-df[l]);
      REAL(EF)[l] += pchisq(n*pow(lambda[l]*gm[j],2.0)*alpha*m[j-1]/tauSq, pow(gm[j], 2.0), 0, 0);
    }
  }
  
  // Return
  UNPROTECT(1);
  return(EF);
}