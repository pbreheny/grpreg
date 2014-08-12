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
SEXP cleanupG(double *a, double *r, int *e, SEXP beta, SEXP iter, SEXP df, SEXP loss);

// Group descent update
void gd_gaussian(double *b, double *x, double *r, int g, int *K1, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a) {
  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + a[j];
  double z_norm = norm(z,K);

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(z_norm, lam1) / (1+lam2);
  if (strcmp(penalty, "grMCP")==0) len = F(z_norm, lam1, lam2, gamma);
  if (strcmp(penalty, "grSCAD")==0) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 | a[K1[g]] != 0) {
    // If necessary, update beta and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      b[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = b[l*p+j]-a[j];
      for (int i=0; i<n; i++) r[i] -= x[n*j+i] * shift;
    }
  }

  // Update df
  if (len > 0) REAL(df)[l] += K * len / z_norm;
  Free(z);
}

SEXP gdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP user_) {

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
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta, iter, df, loss;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  PROTECT(loss = allocVector(REALSXP, L));

  // Intermediate quantities
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int converged, lstart, ng, nv, violations;
  double shift, l1, l2;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    lstart = 1;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    if (l != 0) {
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
	res = cleanupG(a, r, e, beta, iter, df, loss);
	return(res);
      }
    }

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
	converged = 0;
	INTEGER(iter)[l]++;
	REAL(df)[l] = 0;

	// Update unpenalized covariates
	for (int j=0; j<K0; j++) {
	  shift = crossprod(X, r, n, j)/n;
	  b[l*p+j] = shift + a[j];
	  for (int i=0; i<n; i++) r[i] -= shift * X[n*j+i];
	  REAL(df)[l] += 1;
	}

	// Update penalized groups
	for (int g=0; g<J; g++) {
	  l1 = lam[l] * m[g] * alpha;
	  l2 = lam[l] * m[g] * (1-alpha);
	  if (e[g]) gd_gaussian(b, X, r, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
	}

	// Check convergence
	if (checkConvergence(b, a, eps, l, p)) {
	  converged  = 1;
	  REAL(loss)[l] = gLoss(r,n);
	  break;
	}
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
	if (e[g]==0) {
	  l1 = lam[l] * m[g] * alpha;
	  l2 = lam[l] * m[g] * (1-alpha);
	  gd_gaussian(b, X, r, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
	  if (b[l*p+K1[g]] != 0) {
	    e[g] = 1;
	    violations++;
	  }
	}
      }

      if (violations==0) {
	REAL(loss)[l] = gLoss(r, n);
	break;
      }
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  res = cleanupG(a, r, e, beta, iter, df, loss);
  return(res);
}
