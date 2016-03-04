#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double crossprod(double *x, double *y, int n, int j);
double sum(double *x, int n);
double max(double *x, int n);
double norm(double *x, int p);
double S(double z, double l);
double F(double z, double l1, double l2, double gamma);
double Fs(double z, double l1, double l2, double gamma);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);
SEXP cleanupB(double *a, double *r, int *e, double *eta, SEXP beta0, SEXP beta, SEXP iter, SEXP df, SEXP Dev);

// Group descent update -- poisson
void gd_poisson(double *b, double *x, double *r, double v, double *eta, int g, int *K1, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a) {

  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + a[j];
  double z_norm = norm(z,K);

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(v * z_norm, lam1) / (v * (1 + lam2));
  if (strcmp(penalty, "grMCP")==0) len = F(v * z_norm, lam1, lam2, gamma) / v;
  if (strcmp(penalty, "grSCAD")==0) len = Fs(v * z_norm, lam1, lam2, gamma) / v;
  if (len != 0 | a[K1[g]] != 0) {
    // If necessary, update b and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      b[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = b[l*p+j]-a[j];
      for (int i=0; i<n; i++) {
	double si = shift*x[j*n+i];
	r[i] -= si;
	eta[i] += si;
      }
    }
  }

  // Update df
  if (len > 0) REAL(df)[l] += K * len / z_norm;
  Free(z);
}

SEXP gdfit_poisson(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_) {

  // Lengths/dimensions
  int n = length(y_);
  int L = length(lambda);
  int J = length(K1_) - 1;
  int p = length(X_)/n;

  // Pointers
  double *X = REAL(X_);
  double *y = REAL(y_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  int *K1 = INTEGER(K1_);
  int K0 = INTEGER(K0_)[0];
  double *lam = REAL(lambda);
  double alpha = REAL(alpha_)[0];
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double gamma = REAL(gamma_)[0];
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int warn = INTEGER(warn_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta0, beta, iter, df, Dev;
  PROTECT(beta0 = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(beta0)[i] = 0;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(Dev = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(Dev)[i] = 0;
  double *b0 = REAL(beta0);
  double *b = REAL(beta);

  // Intermediate quantities
  double a0 = 0; // Beta0 from previous iteration
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *eta = Calloc(n, double);
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int converged, lstart, ng, nv, violations;
  double shift, l1, l2, mu, v;

  // Initialization
  double ybar = sum(y, n)/n;
  a0 = b0[0] = log(ybar);
  double nullDev = 0;
  for (int i=0;i<n;i++) if (y[i]!=0) nullDev += y[i]*log(y[i]/ybar);
  for (int i=0; i<n; i++) eta[i] = a0;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Dev)[0] = nullDev;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      a0 = b0[l-1];
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax, gmax
      ng = 0;
      nv = 0;
      for (int g=0; g<J; g++) {
	if (a[K1[g]] != 0) {
	  ng++;
	  nv = nv + (K1[g+1]-K1[g]);
	}
      }
      if (ng > gmax | nv > dfmax) {
	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	res = cleanupB(a, r, e, eta, beta0, beta, iter, df, Dev);
	return(res);
      }
    }

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
	INTEGER(iter)[l]++;
	REAL(Dev)[l] = 0;
	v = exp(max(eta, n));
	for (int i=0; i<n; i++) {
	  mu = exp(eta[i]);
	  r[i] = (y[i] - mu)/v;
	  if (y[i]!=0) REAL(Dev)[l] += y[i]*log(y[i]/mu);
	}

	// Check for saturation
	if (REAL(Dev)[l]/nullDev < .001) {
	  if (warn) warning("Model saturated; exiting...");
	  for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	  res = cleanupB(a, r, e, eta, beta0, beta, iter, df, Dev);
	  return(res);
	}

	// Update intercept
	shift = sum(r, n)/n;
	b0[l] = shift + a0;
	for (int i=0; i<n; i++) {
	  r[i] -= shift;
	  eta[i] += shift;
	}
	REAL(df)[l] = 1;

	// Update unpenalized covariates
	for (int j=0; j<K0; j++) {
	  shift = crossprod(X, r, n, j)/n;
	  b[l*p+j] = shift + a[j];
	  for (int i=0; i<n; i++) {
	    double si = shift * X[n*j+i];
	    r[i] -= si;
	    eta[i] += si;
	  }
	  REAL(df)[l]++;
	}

	// Update penalized groups
	for (int g=0; g<J; g++) {
	  l1 = lam[l] * m[g] * alpha;
	  l2 = lam[l] * m[g] * (1-alpha);
	  if (e[g]) gd_poisson(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
	}

	// Check convergence
	converged = 0;
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
	if (e[g]==0) {
	  l1 = lam[l] * m[g] * alpha;
	  l2 = lam[l] * m[g] * (1-alpha);
	  gd_poisson(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
	  if (b[l*p+K1[g]] != 0) {
	    e[g] = 1;
	    violations++;
	  }
	}
      }

      if (violations==0) break;
      a0 = b0[l];
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  res = cleanupB(a, r, e, eta, beta0, beta, iter, df, Dev);
  return(res);
}
