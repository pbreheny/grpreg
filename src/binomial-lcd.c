#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double crossprod(double *x, double *y, int n, int j);
double sum(double *x, int n);
double norm(double *x, int p);
double S(double z, double l);
double F(double z, double l1, double l2, double gamma);
double Fs(double z, double l1, double l2, double gamma);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);
void cleanupB(double *a, double *r, int *e, double *eta);

// Groupwise local coordinate descent updates -- binomial
void gLCD_binomial(double *b, char *penalty, double *x, double *r, double *eta, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *a, double delta, int *e) {

  // Calculate v
  int K = K1[g+1] - K1[g];
  double v = 0.25;

  // Make initial local approximation
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
      double u = crossprod(x, r, n, j)/n + a[j];
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau*v/lam1*sG);
	if (strcmp(penalty, "gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
      }
      b[l*p+j] = S(v*u, ljk) / (v*(1+lam2));

      // Update r, eta, sG, df
      double shift = b[l*p+j] - a[j];
      if (shift != 0) {
	for (int i=0; i<n; i++) {
	  double si = shift*x[j*n+i];
	  r[i] -= si;
	  eta[i] += si;
	}
	if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }
      df[l] = df[l] + fabs(b[l*p+j]) / fabs(u);
    }
  }
}

// KKT check
int gLCD_bCheck(double *b, char *penalty, double *x, double *r, double *eta, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *a,  int *e) {

  // Make initial local approximation
  int violations = 0;
  int K = K1[g+1] - K1[g];
  double v = 0.25;
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
      double u = crossprod(x, r, n, j)/n + a[j];
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau*v/lam1*sG);
      }

      // Update if necessary
      if (v*fabs(u) > ljk) {
	e[j] = 1;
	violations++;
	b[l*p+j] = S(v*u, ljk) / (v*(1+lam2));
	for (int i=0; i<n; i++) {
	  double si = b[l*p+j] * x[j*n+i];
	  r[i] -= si;
	  eta[i] += si;
	}
	if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }
    }
  }
  return(violations);
}

void gpPathFit_binomial(double *b0, double *b, int *iter, double *df, double *Dev, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, double *group_multiplier, int *dfmax_, int *gmax_, int *warn_, int *user_)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double delta=delta_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; int gmax=gmax_[0]; int warn = warn_[0]; int user = user_[0];
  double a0=0;
  double *r = Calloc(n, double);
  double *eta = Calloc(n, double);
  /* double *w = Calloc(n, double); */
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double shift, si;
  int lstart, violations;

  // Initialization
  double ybar = sum(y,n)/n;
  a0 = b0[0] = log(ybar/(1-ybar));
  double nullDev = 0;
  for (int i=0;i<n;i++) nullDev += - y[i]*log(ybar) - (1-y[i])*log(1-ybar);
  for (int i=0; i<n; i++) eta[i] = a0;
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=0; j<p; j++) {
      double z=0;
      for (int i=0; i<n; i++) z += 0.25 * x[j*n+i] * (y[i]-a0);
      a[j] = z/n;
      e[j] = 1;
      for (int i=0; i<n; i++) {
	si = a[j] * x[j*n+i];
	r[i] -= si;
	eta[i] += si;
      }
    }
  }

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user | strcmp(penalty, "gBridge")==0) {
    lstart = 0;
  } else {
    lstart = 1;
    Dev[0] = nullDev;
  }

  // Path
  double pi;
  for (int l=lstart; l<L; l++) {
    if (l != 0) {
      a0 = b0[l-1];
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
	cleanupB(a, r, e, eta);
	return;
      }
    }

    while (iter[l] < max_iter) {
      while (iter[l] < max_iter) {
	int converged = 0;
	iter[l]++;

	// Approximate L
	Dev[l] = 0;
	for (int i=0; i<n; i++) {
	  if (eta[i] > 10) {
	    pi = 1;
	  } else if (eta[i] < -10) {
	    pi = 0;
	  } else {
	    pi = exp(eta[i])/(1+exp(eta[i]));
	  }
	  r[i] = (y[i] - pi) / 0.25;
	  Dev[l] += - y[i]*log(pi) - (1-y[i])*log(1-pi);
	}

	// Check for saturation
	if (Dev[l]/nullDev < .01) {
	  if (warn) warning("Model saturated; exiting...");
	  for (int ll=l; ll<L; ll++) iter[ll] = NA_INTEGER;
	  cleanupB(a, r, e, eta);
	  return;
	}

	// Update intercept
	shift = sum(r, n)/n;
	b0[l] = shift + a0;
	for (int i=0; i<n; i++) {
	  r[i] -= shift;
	  eta[i] += shift;
	}
	df[l] = 1;
  
	// Update unpenalized covariates
	for (int j=0; j<K0; j++) {
	  shift = crossprod(x, r, n, j)/n;
	  b[l*p+j] = shift + a[j];
	  for (int i=0; i<n; i++) {
	    double si = shift * x[n*j+i];
	    r[i] -= si;
	    eta[i] += si;
	  }
	  df[l] = df[l] + 1;
	}

	// Update penalized groups
	for (int g=0; g<J; g++) {
	  gLCD_binomial(b, penalty, x, r, eta, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, df, a, delta, e);
	}

	// Check convergence
	if (checkConvergence(b, a, eps, l, p)) {
	  converged  = 1;
	  break;
	}
	a0 = b0[l];
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
	violations += gLCD_bCheck(b, penalty, x, r, eta, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, a, e);
      }

      if (violations==0) break;
      a0 = b0[l];
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  cleanupB(a, r, e, eta);
  return;
}
