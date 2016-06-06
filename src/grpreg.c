#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
SEXP gdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP user_);
SEXP gdfit_binomial(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP gdfit_poisson(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP gdfit_cox(SEXP X_, SEXP d_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP lcdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP user_);
SEXP lcdfit_binomial(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP lcdfit_poisson(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP lcdfit_cox(SEXP X_, SEXP d_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_);
SEXP standardize(SEXP X_);
SEXP maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_);

// Cleanup
SEXP cleanupG(double *a, double *r, int *e, SEXP beta, SEXP iter, SEXP df, SEXP loss) {
  Free(a);
  Free(r);
  Free(e);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, iter);
  SET_VECTOR_ELT(res, 2, df);
  SET_VECTOR_ELT(res, 3, loss);
  UNPROTECT(5);
  return(res);
}
SEXP cleanupB(double *a, double *r, int *e, double *eta, SEXP beta0, SEXP beta, SEXP iter, SEXP df, SEXP Dev) {
  Free(a);
  Free(r);
  Free(e);
  Free(eta);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, df);
  SET_VECTOR_ELT(res, 4, Dev);
  UNPROTECT(6);
  return(res);
}
// Memory handling and output formatting, Cox
SEXP cleanupCox(double *h, double *a, double *r, int *e, double *eta, double *haz, double *rsk, SEXP beta, SEXP Dev, SEXP iter, SEXP Eta, SEXP df) {
  Free(h);
  Free(a);
  Free(r);
  Free(e);
  Free(eta);
  Free(haz);
  Free(rsk);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, iter);
  SET_VECTOR_ELT(res, 2, df);
  SET_VECTOR_ELT(res, 3, Dev);
  SET_VECTOR_ELT(res, 4, Eta);
  UNPROTECT(6);
  return(res);
}

// Check for convergence of beta[l]
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J) {
  int j;
  int converged = 1;
  for (j=0; j < J; j++) {
    if (beta[l*J+j]!=0 & beta_old[j]!=0) {
      if (fabs(beta[l*J+j]-beta_old[j]) > eps) {
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

static R_CallMethodDef callMethods[] = {
  {"gdfit_gaussian", (DL_FUNC) &gdfit_gaussian, 14},
  {"gdfit_binomial", (DL_FUNC) &gdfit_binomial, 15},
  {"gdfit_poisson", (DL_FUNC) &gdfit_poisson, 15},
  {"gdfit_cox", (DL_FUNC) &gdfit_cox, 15},
  {"lcdfit_gaussian", (DL_FUNC) &lcdfit_gaussian, 16},
  {"lcdfit_binomial", (DL_FUNC) &lcdfit_binomial, 17},
  {"lcdfit_poisson", (DL_FUNC) &lcdfit_poisson, 17},
  {"lcdfit_cox", (DL_FUNC) &lcdfit_cox, 17},
  {"standardize", (DL_FUNC) &standardize, 1},
  {"maxprod", (DL_FUNC) &maxprod, 4},
  {NULL, NULL, 0}
};

void R_init_grpreg(DllInfo *info) {
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
