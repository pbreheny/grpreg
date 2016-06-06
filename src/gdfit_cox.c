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
SEXP cleanupCox(double *h, double *a, double *r, int *e, double *eta, double *haz,
		double *rsk, SEXP beta, SEXP Dev, SEXP iter, SEXP Eta, SEXP df);

// Group descent update -- cox
void gd_cox(double *b, double *x, double *r, double *eta, double v, int g,
	    int *K1, int n, int l, int p, const char *penalty, double lam1,
	    double lam2, double gamma, SEXP df, double *a) {

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

SEXP gdfit_cox(SEXP X_, SEXP d_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_) {

  // Lengths/dimensions
  int n = length(d_);
  int L = length(lambda);
  int J = length(K1_) - 1;
  int p = length(X_)/n;

  // Pointers
  double *X = REAL(X_);
  double *d = REAL(d_);
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
  SEXP res, beta, Loss, iter, df, Eta;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(Loss = allocVector(REALSXP, L));
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int i=0; i<(L*n); i++) REAL(Eta)[i] = 0;
  double *b = REAL(beta);
  double *ETA = REAL(Eta);

  // Intermediate quantities
  double *a = Calloc(p, double);  // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double *r = Calloc(n, double);
  double *h = Calloc(n, double);  
  double *haz = Calloc(n, double);
  double *rsk = Calloc(n, double);  
  double *eta = Calloc(n, double);
  int *e = Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int converged, lstart, ng, nv, violations;
  double shift, l1, l2, nullDev, v, s;

  // Initialization
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  rsk[n-1] = 1;
  for (int i=n-2; i>=0; i--) rsk[i] = rsk[i+1] + 1;
  nullDev = 0;
  for (int i=0; i<n; i++) nullDev -= d[i]*log(rsk[i]);
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Loss)[0] = nullDev;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
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
	res = cleanupCox(h, a, r, e, eta, haz, rsk, beta, Loss, iter, Eta, df);
	return(res);
      }
    }

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
	INTEGER(iter)[l]++;
	REAL(Loss)[l] = 0;
	REAL(df)[l] = 0;
	
	// Calculate haz, risk
	for (int i=0; i<n; i++) haz[i] = exp(eta[i]);
	rsk[n-1] = haz[n-1];
	for (int i=n-2; i>=0; i--) {
	  rsk[i] = rsk[i+1] + haz[i];
	}
	for (int i=0; i<n; i++) {
	  REAL(Loss)[l] += d[i]*eta[i] - d[i]*log(rsk[i]);
	}
	  
	// Approximate L
	h[0] = d[0]/rsk[0];
	v = 1;
	for (int i=1; i<n; i++) {
	  h[i] = h[i-1] + d[i]/rsk[i];
	}
	for (int i=0; i<n; i++) {
	  h[i] = h[i]*haz[i];
	  s = d[i] - h[i];
	  if (h[i]==0) r[i]=0;
	  else r[i] = s/v;
	}

	// Check for saturation
	if (REAL(Loss)[l]/nullDev < .01) {
	  if (warn) warning("Model saturated; exiting...");
	  for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	  res = cleanupCox(h, a, r, e, eta, haz, rsk, beta, Loss, iter, Eta, df);
	  return(res);
	}

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
	  if (e[g]) gd_cox(b, X, r, eta, v, g, K1, n, l, p, penalty, l1, l2,
			   gamma, df, a);
	}

	// Check convergence
	converged = checkConvergence(b, a, eps, l, p);
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
	if (converged) break;
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
	if (e[g]==0) {
	  l1 = lam[l] * m[g] * alpha;
	  l2 = lam[l] * m[g] * (1-alpha);
	  gd_cox(b, X, r, eta, v, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
	  if (b[l*p+K1[g]] != 0) {
	    e[g] = 1;
	    violations++;
	  }
	}
      }

      if (violations==0) {
	for (int i=0; i<n; i++) ETA[l*n+i] = eta[i];
	break;
      }
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  res = cleanupCox(h, a, r, e, eta, haz, rsk, beta, Loss, iter, Eta, df);
  return(res);
}
