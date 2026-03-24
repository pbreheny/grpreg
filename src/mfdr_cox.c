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

SEXP mfdr_cox (SEXP fit) {
  
  // Declarations
  int n = INTEGER(getListElement(fit, "n"))[0];
  int L = ncols(getListElement(fit, "beta"));
  int ng = length(getListElement(fit, "group.multiplier"));  
  int *gl= INTEGER(getListElement(fit, "gl"));  
  int ck;
  double *Eta = REAL(getListElement(fit, "linear.predictors")); 
  double *d = REAL(getListElement(fit, "fail"));
  double *X = REAL(getListElement(fit, "XX"));
  double *gm = REAL(getListElement(fit, "group.multiplier"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double tauSq;
  double *w = R_Calloc(n, double);
  double *haz = R_Calloc(n, double);
  double *rsk = R_Calloc(n, double);
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;
  
  // Calculation
  for (int l=0; l<L; l++) {
    
    //Calculate Risk
    for (int i=0; i<n; i++){
      haz[i] = exp(Eta[n*l + i]);
    }
      rsk[n-1] = haz[n-1];
    for (int i=n-2; i>=0; i--){
      rsk[i] = rsk[i+1] + haz[i];  
    }
    
    // Find W for current lam
    for (int j=0; j<n; j++){
      w[j] = 0;
      for (int i=0; i<=j; i++){
        w[j] += d[i]*haz[j]/rsk[i]*(1-haz[j]/rsk[i]);
        }
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
  free(w);
  free(haz);
  free(rsk);
  UNPROTECT(1);
  return(EF);
}