#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
void gpPathFit_gaussian(double *beta, int *iter, double *df, double *loss, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, int *dfmax_, int *gmax_, double *group_multiplier, int *user_);
void gpPathFit_binomial(double *beta0, double *beta, int *iter, double *df, double *Dev, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, double *group_multiplier, int *dfmax_, int *gmax_, int *warn_, int *user_);
void grPathFit_gaussian(double *beta, int *iter, double *df, double *loss, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *group_multiplier, int *dfmax_, int *gmax_, int *user_);
void grPathFit_binomial(double *beta0, double *beta, int *iter, double *df, double *Dev, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *group_multiplier, int *dfmax_, int *gmax_, int *warn_, int *user_);
SEXP standardize(SEXP X_);
SEXP maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_);

// Cleanup
void cleanupG(double *a, double *r, int *e) {
  Free(a);
  Free(r);
  Free(e);
}
void cleanupB(double *a, double *r, int *e, double *eta) {
  Free(a);
  Free(r);
  Free(e);
  Free(eta);
}

// Check for convergence of beta[l]
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J) {
  int j;
  int converged = 1;
  for (j=0; j < J; j++) {
    if (beta[l*J+j]!=0 & beta_old[j]!=0) {
      if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps) {
	converged = 0;
	break;
      }
    } else if (beta[l*J+j]==0 & beta_old[j]!=0) {
      converged = 0;
      break;
    } else if (beta[l*J+j]!=0 & beta_old[j]==0) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

// Cross product of the jth column of x with y
double crossprod(double *x, double *y, int n, int j) {
  double val = 0;
  int nn = n*j;
  for (int i=0; i<n; i++) val += x[nn+i] * y[i];
  return(val);
}

// Sum of x
double sum(double *x, int n) {
  double val = 0;
  for (int i=0; i<n; i++) val += x[i];
  return(val);
}

// Gaussian loss
double gLoss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
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

static R_CMethodDef cMethods[] = {
  {"gpPathFit_gaussian", (DL_FUNC) &gpPathFit_gaussian, 24},
  {"gpPathFit_binomial", (DL_FUNC) &gpPathFit_binomial, 26},
  {"grPathFit_gaussian", (DL_FUNC) &grPathFit_gaussian, 22},
  {"grPathFit_binomial", (DL_FUNC) &grPathFit_binomial, 24},
  {NULL, NULL, 0}
};

static R_CallMethodDef callMethods[] = {
  {"standardize", (DL_FUNC) &standardize, 1},
  {"maxprod", (DL_FUNC) &maxprod, 4},
  {NULL, NULL, 0}
};

void R_init_grpreg(DllInfo *info) {
  R_registerRoutines(info,cMethods,callMethods,NULL,NULL);
}
