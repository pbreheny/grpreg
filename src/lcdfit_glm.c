#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *x, double *y, int n, int j);
double sum(double *x, int n);
double max(double *x, int n);
double p_binomial(double eta);
double norm(double *x, int p);
double S(double z, double l);
double F(double z, double l1, double l2, double gamma);
double Fs(double z, double l1, double l2, double gamma);
double MCP(double theta, double l, double a);
double dMCP(double theta, double l, double a);

// Groupwise local coordinate descent updates
void gLCD_glm(double *b, const char *penalty, double *x, double *r, double v, double *eta, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, SEXP df, double *a, double delta, int *e, double *maxChange) {

  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  double shift;
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "gel")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j])/v;
  if (strcmp(penalty, "cMCP")==0) {
    lam1 = sqrt(lam1);
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(a[j]/v, lam1, gamma);
  }
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(a[j]/v);
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
      double u = crossprod(x, r, n, j)/n + a[j];
      double ljk=0;
      if (lam1 != 0) {
        if (strcmp(penalty, "cMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(b[l*p+j], lam1, gamma);
        if (strcmp(penalty, "gel")==0) ljk = lam1*exp(-tau/lam1*sG);
        if (strcmp(penalty, "gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
      }
      b[l*p+j] = S(v*u, ljk) / (v*(1+lam2));

      // Update r, eta, sG, df
      shift = b[l*p+j] - a[j];
      if (shift != 0) {
        if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
        for (int i=0; i<n; i++) {
          double si = shift*x[j*n+i];
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
}

// KKT check
int gLCD_glm_check(double *b, const char *penalty, double *x, double *r, double v, double *eta, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *a,  int *e) {

  // Make initial local approximation
  int violations = 0;
  int K = K1[g+1] - K1[g];
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
      double u = crossprod(x, r, n, j)/n;
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

SEXP lcdfit_glm(SEXP X_, SEXP y_, SEXP family_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP delta_, SEXP gamma_, SEXP tau_, SEXP max_iter_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_) {

  // Lengths/dimensions
  int n = length(y_);
  int L = length(lambda);
  int J = length(K1_) - 1;
  int p = length(X_)/n;

  // Pointers
  double *X = REAL(X_);
  double *y = REAL(y_);
  const char *family = CHAR(STRING_ELT(family_, 0));
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
  int warn = INTEGER(warn_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta0, beta, iter, df, Dev;
  PROTECT(res = allocVector(VECSXP, 5));
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
  double a0;   // Beta0 from previous iteration
  double *a = Calloc(p, double);
  double *r = Calloc(n, double);
  int *e = Calloc(p, int);
  double *eta = Calloc(n, double);
  int lstart, ng, nv, violations;
  double shift, l1, l2, mu, v, maxChange;

  // Initialization
  double ybar = sum(y,n)/n;
  double nullDev = 0;
  if (strcmp(family, "binomial")==0) {
    a0 = b0[0] = log(ybar/(1-ybar));
    for (int i=0; i<n; i++) nullDev -= 2*y[i]*log(ybar) + 2*(1-y[i])*log(1-ybar);
  } else if (strcmp(family, "poisson")==0) {
    a0 = b0[0] = log(ybar);
    for (int i=0;i<n;i++) {
      if (y[i]!=0) nullDev += 2*(y[i]*log(y[i]/ybar) + ybar - y[i]);
      else nullDev += 2*ybar;
    }
  }
  for (int i=0; i<n; i++) eta[i] = a0;
  if (strcmp(penalty, "gBridge")==0) {
    for (int i=0; i<n; i++) r[i] = y[i] - a0;
    for (int j=0; j<p; j++) {
      double z=0;
      for (int i=0; i<n; i++) z += 0.25 * X[j*n+i] * (y[i] - a0);
      a[j] = z/n;
      e[j] = 1;
      for (int i=0; i<n; i++) {
        double si = a[j] * X[j*n+i];
        r[i] -= si;
        eta[i] += si;
      }
    }
  } else {
    for (int j=0; j<p; j++) a[j] = 0;
    for (int j=0; j<p; j++) e[j] = 0;
  }

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

        // Approximate L
        REAL(Dev)[l] = 0;
        if (strcmp(family, "binomial")==0) {
          v = 0.25;
          for (int i=0; i<n; i++) {
            mu = p_binomial(eta[i]);
            r[i] = (y[i] - mu) / v;
            if (y[i]==1) REAL(Dev)[l] -= 2*log(mu);
            else REAL(Dev)[l] -= 2*log(1-mu);
          }
        } else if (strcmp(family, "poisson")==0) {
          v = exp(max(eta, n));
          for (int i=0; i<n; i++) {
            mu = exp(eta[i]);
            r[i] = (y[i] - mu)/v;
            if (y[i]!=0) REAL(Dev)[l] += 2*(y[i]*log(y[i]/mu) + mu - y[i]);
            else REAL(Dev)[l] += 2*mu;
          }
        }

        // Check for saturation
        if (REAL(Dev)[l]/nullDev < .01) {
          if (warn) warning("Model saturated; exiting...");
          for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }

        // Update intercept
        shift = sum(r, n)/n;
        b0[l] = shift + a0;
        for (int i=0; i<n; i++) {
          r[i] -= shift;
          eta[i] += shift;
        }
        REAL(df)[l] = 1;
        maxChange = fabs(shift);

        // Update unpenalized covariates
        for (int j=0; j<K0; j++) {
          shift = crossprod(X, r, n, j)/n;
          if (fabs(shift) > maxChange) maxChange = fabs(shift);
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
          gLCD_glm(b, penalty, X, r, v, eta, g, K1, n, l, p, l1, l2, gamma, tau, df, a, delta, e, &maxChange);
        }

        // Check convergence
        a0 = b0[l];
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
        l1 = lam[l] * m[g] * alpha;
        l2 = lam[l] * m[g] * (1-alpha);
        violations += gLCD_glm_check(b, penalty, X, r, v, eta, g, K1, n, l, p, l1, l2, gamma, tau, a, e);
      }

      if (violations==0) break;
      a0 = b0[l];
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  Free(a);
  Free(r);
  Free(e);
  Free(eta);
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, df);
  SET_VECTOR_ELT(res, 4, Dev);
  UNPROTECT(6);
  return(res);
}
