#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP gdfit_glm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gdfit_cox(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gdfit_gaussian(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcdfit_glm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcdfit_cox(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lcdfit_gaussian(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP maxgrad(SEXP, SEXP, SEXP, SEXP);
extern SEXP maxprod(SEXP, SEXP, SEXP, SEXP);
extern SEXP standardize(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"gdfit_glm",                (DL_FUNC) &gdfit_glm,                16},
  {"gdfit_cox",                (DL_FUNC) &gdfit_cox,                15},
  {"gdfit_gaussian",           (DL_FUNC) &gdfit_gaussian,           15},
  {"lcdfit_glm",               (DL_FUNC) &lcdfit_glm,               18},
  {"lcdfit_cox",               (DL_FUNC) &lcdfit_cox,               17},
  {"lcdfit_gaussian",          (DL_FUNC) &lcdfit_gaussian,          16},
  {"maxgrad",                  (DL_FUNC) &maxgrad,                   4},
  {"maxprod",                  (DL_FUNC) &maxprod,                   4},
  {"standardize",              (DL_FUNC) &standardize,               1},
  {NULL, NULL, 0}
};

void R_init_grpreg(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

// Cross product of the jth column of x with y
double crossprod(double *x, double *y, int n, int j) {
  double val = 0;
  int nn = n*j;
  for (int i=0; i<n; i++) val += x[nn+i] * y[i];
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of x
double sum(double *x, int n) {
  double val = 0;
  for (int i=0; i<n; i++) val += x[i];
  return(val);
}

// Max of x
double max(double *x, int n) {
  double val = x[0];
  for (int i=1; i<n; i++) {
    if (x[i] > val) val = x[i];
  }
  return(val);
}

// Gaussian loss
double gLoss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Pr(y=1) for binomial
double p_binomial(double eta) {
  if (eta > 10) {
    return(1);
  } else if (eta < -10) {
    return(0);
  } else {
    return(exp(eta)/(1+exp(eta)));
  }
}

// Euclidean norm
double norm(double *x, int p) {
  double x_norm = 0;
  for (int j=0; j<p; j++) x_norm = x_norm + pow(x[j],2);
  x_norm = sqrt(x_norm);
  return(x_norm);
}

// Soft-thresholding operator
double S(double z, double l) {
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

// Firm-thresholding operator
double F(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(1+l2-1/gamma));
  else return(z/(1+l2));
}

// SCAD-modified firm-thresholding operator
double Fs(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(1+l2));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(1-1/(gamma-1)+l2));
  else return(z/(1+l2));
}

// MCP penalty
double MCP(double theta, double l, double a) {
  theta = fabs(theta);
  if (theta <= a*l) return(l*theta - pow(theta,2)/(2*a));
  else return(a*pow(l,2)/2);
}

// MCP penalization rate
double dMCP(double theta, double l, double a) {
  theta = fabs(theta);
  if (theta < a*l) return(l-theta/a);
  else return(0);
}
