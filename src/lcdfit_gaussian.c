#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *x, double *y, int n, int j);
double norm(double *x, int p);
double S(double z, double l);
double F(double z, double l1, double l2, double gamma);
double Fs(double z, double l1, double l2, double gamma);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);
double gLoss(double *r, int n);

// Groupwise local coordinate descent updates
void gLCD_gaussian(double *b, const char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, SEXP df, double *a, double delta, int *e, double *maxChange) {

  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  double sG = 0; // Sum of inner penalties for group
  double shift;
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
        shift = b[l*p+j] - a[j];
        if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
        for (int i=0; i<n; i++) r[i] = r[i] - shift * x[n*j+i];
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
      shift = b[l*p+j] - a[j];
      if (shift != 0) {
        if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
        for (int i=0; i<n; i++) r[i] -= shift*x[n*j+i];
        if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
        if (strcmp(penalty, "gel")==0) sG = sG + fabs(b[l*p+j]) - fabs(a[j]);
        if (strcmp(penalty, "cMCP")==0) sG = sG + MCP(b[l*p+j], lam1, gamma) - MCP(a[j], lam1, gamma);
      }

      // Update df
      REAL(df)[l] += fabs(b[l*p+j]) / fabs(z);
    }
  }
}

// Check inactive variables
int gLCD_gCheck(double *b, const char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *a, int *e) {

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

SEXP lcdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP user_) {

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
  double delta = REAL(delta_)[0];
  double gamma = REAL(gamma_)[0];
  double tau = REAL(tau_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta, loss, Eta, df, iter;
  PROTECT(res = allocVector(VECSXP, 5));
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(loss = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(loss)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;

  // Intermediate quantities
  double *a = Calloc(p, double);
  double *r = Calloc(n, double);
  int *e = Calloc(p, int);
  if (strcmp(penalty, "gBridge")==0) {
    for (int i=0; i<n; i++) r[i] = y[i];
    for (int j=0; j<p; j++) {
      a[j] = crossprod(X, r, n, j)/n;
      for (int i=0; i<n; i++) r[i] = r[i] - a[j] * X[j*n+i];
      e[j] = 1;
    }
  } else {
    for (int j=0; j<p; j++) a[j] = 0;
    for (int j=0; j<p; j++) e[j] = 0;
    for (int i=0; i<n; i++) r[i] = y[i];
  }
  int lstart, ng, nv, violations;
  double shift, l1, l2, maxChange;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  double rss = gLoss(r,n);
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(loss)[0] = rss;
  }
  double sdy = sqrt(rss/n);

  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
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
      if (ng > gmax || nv > dfmax || tot_iter == max_iter) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }
    }

    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        INTEGER(iter)[l]++;
        tot_iter++;
        REAL(df)[l] = 0;

        // Update unpenalized covariates
        maxChange = 0;
        for (int j=0; j<K0; j++) {
          shift = crossprod(X, r, n, j)/n;
          if (fabs(shift) > maxChange) maxChange = fabs(shift);
          b[l*p+j] = shift + a[j];
          for (int i=0; i<n; i++) r[i] -= shift * X[n*j+i];
          REAL(df)[l]++;
        }

        // Update penalized groups
        for (int g=0; g<J; g++) {
          l1 = lam[l] * m[g] * alpha;
          l2 = lam[l] * m[g] * (1-alpha);
          gLCD_gaussian(b, penalty, X, r, g, K1, n, l, p, l1, l2, gamma, tau, df, a, delta, e, &maxChange);
        }

        // Check for convergence      
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps*sdy) break;
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
        l1 = lam[l] * m[g] * alpha;
        l2 = lam[l] * m[g] * (1-alpha);
        violations += gLCD_gCheck(b, penalty, X, r, g, K1, n, l, p, l1, l2, gamma, tau, a, e);
      }

      if (violations==0) {
        REAL(loss)[l] = gLoss(r, n);
        for (int i=0; i<n; i++) REAL(Eta)[n*l+i] = y[i] - r[i];
        break;
      }
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  Free(a);
  Free(r);
  Free(e);
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, Eta);
  SET_VECTOR_ELT(res, 3, df);
  SET_VECTOR_ELT(res, 4, iter);
  UNPROTECT(6);
  return(res);
}
