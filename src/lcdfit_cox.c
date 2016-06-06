#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double crossprod(double *x, double *y, int n, int j);
double wcrossprod(double *X, double *y, double *w, int n, int j);
double wsqsum(double *X, double *w, int n, int j);
double norm(double *x, int p);
double S(double z, double l);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);
SEXP cleanupCox(double *h, double *a, double *r, int *e, double *eta, double *haz,
		double *rsk, SEXP beta, SEXP Dev, SEXP iter, SEXP Eta, SEXP df);

// Groupwise local coordinate descent updates -- Cox
void gLCD_cox(double *b, const char *penalty, double *X, double *r, double *eta,
	      double *h, int g, int *K1, int n, int l, int p, double lam1, double lam2,
	      double gamma, double tau, SEXP df, double *a, double delta, int *e) {

  // Pre-calculcate v
  int K = K1[g+1] - K1[g];
  double xwr, u;
  double *v = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]) {
      v[j-K1[g]] = wsqsum(X, h, n, j)/n;
    } else {
      v[j-K1[g]] = 1;
    }
  }

  // Make initial local approximation
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "gel")==0) {
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j])/v[j-K1[g]];
  }
  if (strcmp(penalty, "cMCP")==0) {
    lam1 = sqrt(lam1);
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(a[j]/v[j-K1[g]], lam1, gamma);
  }
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j])/v[j-K1[g]];
    if (sG==0) return;
    if (sG < delta) {
      for (int j=K1[g]; j<K1[g+1]; j++) {
	b[l*p+j] = 0;
	for (int i=0; i<n; i++) r[i] = r[i] - (b[l*p+j] - a[j]) * X[n*j+i];
      }
      return;
    }
  }

  // Coordinate descent
  for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]) {

      // Calculate u, v
      xwr = wcrossprod(X, r, h, n, j);
      u = xwr/n + v[j-K1[g]]*a[j];

      // Update b
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau/lam1*sG);
	if (strcmp(penalty, "gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
      }
      b[l*p+j] = S(u, ljk) / (v[j-K1[g]]*(1+lam2));

      // Update r, eta, sG, df
      double shift = b[l*p+j] - a[j];
      if (shift != 0) {
	for (int i=0; i<n; i++) {
	  double si = shift*X[j*n+i];
	  r[i] -= si;
	  eta[i] += si;
	}
	if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
	if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }
      REAL(df)[l] += fabs(b[l*p+j]) / fabs(u);
    }
  }
  Free(v);
}

// KKT check
int gLCD_cCheck(double *b, const char *penalty, double *X, double *r, double *eta,
		double *h, int g, int *K1, int n, int l, int p, double lam1, double lam2,
		double gamma, double tau, double *a,  int *e) {

  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  double *v = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]) {
      v[j-K1[g]] = wsqsum(X, h, n, j)/n;
    } else {
      v[j-K1[g]] = 1;
    }
  }
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "gel")==0) {
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j])/v[j-K1[g]];
  }
  if (strcmp(penalty, "cMCP")==0) {
    lam1 = sqrt(lam1);
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(a[j]/v[j-K1[g]], lam1, gamma);
  }
  Free(v);

  // Check
 int violations = 0;
 for (int j=K1[g]; j<K1[g+1]; j++) {
    if (e[j]==0) {
      double xwr = wcrossprod(X, r, h, n, j)/n;
      double ljk=0;
      if (lam1 != 0) {
	if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
	if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau/lam1*sG);
      }
      if (fabs(xwr) > ljk) {
	e[j] = 1;
	violations++;
      }
    }
  }
  return(violations);
}

SEXP lcdfit_cox(SEXP X_, SEXP d_, SEXP penalty_, SEXP K1_, SEXP K0_,
		SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_,
		SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_,
		SEXP gmax_, SEXP warn_, SEXP user_) {

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
  double delta = REAL(delta_)[0];
  double gamma = REAL(gamma_)[0];
  double tau = REAL(tau_)[0];
  int max_iter = INTEGER(max_iter_)[0];
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
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  double *r = Calloc(n, double);
  double *h = Calloc(n, double);
  double *haz = Calloc(n, double);
  double *rsk = Calloc(n, double);  
  double *eta = Calloc(n, double);
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  int converged, lstart, ng, nv, violations;
  double shift, l1, l2, nullDev, u, v, s, xwr, xwx;

  // Initialization
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
  if (strcmp(penalty, "gBridge")==0) {
    h[0] = d[0]/rsk[0];
    for (int i=1; i<n; i++) {
      h[i] = h[i-1] + d[i]/rsk[i];
    }
    for (int i=0; i<n; i++) {
      s = d[i] - h[i];
      if (h[i]==0) r[i]=0;
      else r[i] = s;
    }
    for (int j=0; j<p; j++) {
      a[j] = crossprod(X, r, n, j)/n;
      e[j] = 1;
      for (int i=0; i<n; i++) eta[i] += a[j] * X[j*n+i];
    }
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
	int nv_old = nv;
	for (int j=K1[g]; j<K1[g+1]; j++) {
	  if (a[j] != 0) nv++;
	}
	if (nv != nv_old) ng++;
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
          else r[i] = s/h[i];
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
	  xwr = wcrossprod(X, r, h, n, j);
	  xwx = wsqsum(X, h, n, j);
	  u = xwr/n;
	  v = xwx/n;
	  shift = u/v;
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
	  gLCD_cox(b, penalty, X, r, eta, h, g, K1, n, l, p, l1, l2, gamma, tau, df, a, delta, e);
	}

	// Check convergence
	if (checkConvergence(b, a, eps, l, p)) {
	  converged  = 1;
	  break;
	}
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
	l1 = lam[l] * m[g] * alpha;
	l2 = lam[l] * m[g] * (1-alpha);
	violations += gLCD_cCheck(b, penalty, X, r, eta, h, g, K1, n, l, p, l1, l2, gamma, tau, a, e);
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
