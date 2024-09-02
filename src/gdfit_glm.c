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

// Group descent update
void gd_glm(double *b, double *x, double *r, double v, double *eta, int g, int *K1, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a, double *maxChange) {

  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = R_Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + a[j];
  double z_norm = norm(z,K);

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(v * z_norm, lam1) / (v * (1 + lam2));
  if (strcmp(penalty, "grMCP")==0) len = F(v * z_norm, lam1, lam2, gamma) / v;
  if (strcmp(penalty, "grSCAD")==0) len = Fs(v * z_norm, lam1, lam2, gamma) / v;
  if (len != 0 || a[K1[g]] != 0) {
    // If necessary, update b and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      b[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = b[l*p+j]-a[j];
      if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
      for (int i=0; i<n; i++) {
        double si = shift*x[j*n+i];
        r[i] -= si;
        eta[i] += si;
      }
    }
  }

  // Update df
  if (len > 0) REAL(df)[l] += K * len / z_norm;
  free(z);
}

SEXP gdfit_glm(SEXP X_, SEXP y_, SEXP family_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP warn_, SEXP user_) {

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
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int warn = INTEGER(warn_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta0, beta, Dev, Eta, df, iter;
  PROTECT(res = allocVector(VECSXP, 6));
  PROTECT(beta0 = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(beta0)[i] = 0;
  double *b0 = REAL(beta0);
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(Dev = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(Dev)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;

  // Intermediate quantities
  double a0 = 0; // Beta0 from previous iteration
  double *r = R_Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *eta = R_Calloc(n, double);
  double *a = R_Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = R_Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int lstart, ng, nv, violations;
  double shift, l1, l2, mu, v, maxChange;

  // Initialization
  double ybar = sum(y, n)/n;
  double nullDev = 0;
  if (strcmp(family, "binomial") == 0) {
    a0 = b0[0] = log(ybar/(1-ybar));
    for (int i=0; i<n; i++) nullDev -= 2*y[i]*log(ybar) + 2*(1-y[i])*log(1-ybar);
  } else if (strcmp(family, "poisson") == 0) {
    a0 = b0[0] = log(ybar);
    for (int i=0;i<n;i++) {
      if (y[i]!=0) nullDev += 2*(y[i]*log(y[i]/ybar) + ybar - y[i]);
      else nullDev += 2*ybar;
    }
  }
  for (int i=0; i<n; i++) eta[i] = a0;

  // If lam[0]=lam_max, skip lam[0] -- closed form solution available
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
        if (REAL(Dev)[l]/nullDev < 0.01) {
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
          if (e[g]) gd_glm(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, &maxChange);
        }

        // Check convergence
        a0 = b0[l];
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }

      // Scan for violations
      violations = 0;
      for (int g=0; g<J; g++) {
        if (e[g]==0) {
          l1 = lam[l] * m[g] * alpha;
          l2 = lam[l] * m[g] * (1-alpha);
          gd_glm(b, X, r, v, eta, g, K1, n, l, p, penalty, l1, l2, gamma, df, a, &maxChange);
          if (b[l*p+K1[g]] != 0) {
            e[g] = 1;
            violations++;
          }
        }
      }

      if (violations==0) {
        for (int i=0; i<n; i++) REAL(Eta)[n*l+i] = eta[i];
        break;
      }
      a0 = b0[l];
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  free(a);
  free(r);
  free(e);
  free(eta);
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, Dev);
  SET_VECTOR_ELT(res, 3, Eta);
  SET_VECTOR_ELT(res, 4, df);
  SET_VECTOR_ELT(res, 5, iter);
  UNPROTECT(7);
  return(res);
}
