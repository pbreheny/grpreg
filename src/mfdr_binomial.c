#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double wsqsum(double *X, double *w, int n, int j);
SEXP getListElement(SEXP list, const char *str);
double pchisq(double x, double df, int lower_tail, int give_log);

SEXP mfdr_binomial(SEXP fit) {
  
  // Declarations
  int n = INTEGER(getListElement(fit, "n"))[0];
  int L = ncols(getListElement(fit, "beta"));
  int ng = length(getListElement(fit, "group.multiplier"));  
  int *gl = INTEGER(getListElement(fit, "gl")); 
  int ck;
  double *pi = REAL(getListElement(fit, "P")); 
  double *X = REAL(getListElement(fit, "XX"));
  double *gm = REAL(getListElement(fit, "group.multiplier"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double tauSq;
  double *w = Calloc(n, double);
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;
  
  // Calculation
  for (int l=0; l<L; l++) {

  // Find W for current lam
    for (int i=0; i<n; i++) {
      w[i] = pi[n*l+i]*(1-pi[n*l+i]);
    }
  
  // Find EF contrib of each grp
    ck = 0;
    for (int j=0; j < ng; j++) {
      tauSq = 0;
      for (int k=0; k < gl[j]; k++) {
        tauSq += wsqsum(X, w, n, ck+k)/n;
      }
      REAL(EF)[l] += pchisq(n*pow(lambda[l]*gm[j]*gm[j], 2.0)*alpha/tauSq, gl[j], 0, 0);
      ck++;
    }
  }
  
  // Return
  Free(w);
  UNPROTECT(1);
  return(EF);
}