#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double crossprod(double *x, double *y, int n, int j);
double norm(double *x, int p);
double S(double z, double l);
double F(double z, double l1, double l2, double gamma);
double Fs(double z, double l1, double l2, double gamma);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);
double gLoss(double *r, int n);
double cleanupG(double *a, double *r, int *e);

// Groupwise local coordinate descent updates
void gLCD_gaussian(double *b, char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *a, double delta, int *e)
{
  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "gel")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j]);
  if (strcmp(penalty, "cMCP")==0) {
    lam1 = sqrt(lam1);
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(a[j], lam1, gamma);
  }
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j]);
    if (sG==0) return;
    if (sG < delta) {
      for (int j=K1[g]; j<K1[g+1]; j++) {
	b[l*p+j] = 0;
	for (int i=0; i<n; i++) r[i] = r[i] - (b[l*p+j] - a[j]) * x[n*j+i];
      }
      return;
    }
  }

  // Coordinate descent
  for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]) {

      // Update b
      double z = crossprod(x, r, n, j)/n + a[j];
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau/lam1*sG);
	if (strcmp(penalty, "gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
      }
      b[l*p+j] = S(z, ljk) / (1+lam2);

      // Update r
      double shift = b[l*p+j] - a[j];
      if (shift != 0) {
	for (int i=0; i<n; i++) r[i] -= shift*x[n*j+i];
	if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }

      // Update df
      df[l] = df[l] + fabs(b[l*p+j]) / fabs(z);
    }
  }
}

// Check inactive variables
int gLCD_gCheck(double *b, char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *a, int *e) {

  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  int violations=0;
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "gel")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j]);
  if (strcmp(penalty, "cMCP")==0) {
    lam1 = sqrt(lam1);
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(a[j], lam1, gamma);
  }

  // Check
  for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]==0) {
      // Compare
      double z = crossprod(x, r, n, j)/n;
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau/lam1*sG);
      }
      if (fabs(z) > ljk) {
	e[j] = 1;
	violations++;
	b[l*p+j] = S(z, ljk) / (1+lam2);
	for (int i=0; i<n; i++) r[i] = r[i] - b[l*p+j] * x[n*j+i];
	if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }
    }
  }
  return(violations);
}

void gpPathFit_gaussian(double *b, int *iter, double *df, double *loss, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, int *dfmax_, int *gmax_, double *group_multiplier, int *user_)
{
  // Initialization of variables
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double delta=delta_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; int gmax=gmax_[0]; int user=user_[0];
  int lstart, violations;
  double *r = Calloc(n, double);
  int *e = Calloc(p, int);
  double *a = Calloc(p, double);
  if (strcmp(penalty, "gBridge")==0) {
    for (int i=0; i<n; i++) r[i] = y[i];
    for (int j=0; j<p; j++) {
      a[j] = crossprod(x, r, n, j)/n;
      for (int i=0; i<n; i++) r[i] = r[i] - a[j] * x[j*n+i];
      e[j] = 1;
    }
  } else {
    for (int j=0; j<p; j++) a[j] = 0;
    for (int j=0; j<p; j++) e[j] = 0;
    for (int i=0; i<n; i++) r[i] = y[i];
  }

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user | strcmp(penalty, "gBridge")==0) {
    lstart = 0;
  } else {
    loss[0] = gLoss(r,n);
    lstart = 1;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    if (l != 0) {
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax, gmax
      int ng = 0;
      int nv = 0;
      for (int g=0; g<J; g++) {
	int nv_old = nv;
	for (int j=K1[g]; j<K1[g+1]; j++) {
	  if (a[j] != 0) nv++;
	}
	if (nv != nv_old) ng++;
      }
      if (ng > gmax | nv > dfmax) {
	for (int ll=l; ll<L; ll++) iter[ll] = NA_INTEGER;
	cleanupG(a, r, e);
	return;
      }
    }

    while (iter[l] < max_iter) {
      while (iter[l] < max_iter) {
	int converged = 0;
	iter[l]++;
	df[l] = 0;

	// Update unpenalized covariates
	for (int j=0; j<K0; j++) {
	  double shift = crossprod(x, r, n, j)/n;
	  b[l*p+j] = shift + a[j];
	  for (int i=0; i<n; i++) r[i] -= shift * x[n*j+i];
	  df[l] = df[l] + 1;
	}

	// Update penalized groups
	for (int g=0; g<J; g++) {
	  gLCD_gaussian(b, penalty, x, r, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, df, a, delta, e);
	}

	// Check for convergence      
	if (checkConvergence(b, a, eps, l, p)) {
	  converged  = 1;
	  loss[l] = gLoss(r,n);
	  break;
	}
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
	violations += gLCD_gCheck(b, penalty, x, r, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, a, e);
      }

      if (violations==0) {
	loss[l] = gLoss(r, n);
	break;
      }
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  cleanupG(a, r, e);
}
