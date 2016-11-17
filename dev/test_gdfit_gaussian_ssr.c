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
int sum_rejections(int *x, int n);
void gd_gaussian(double *b, double *x, double *r, int g, int *K1, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a);
SEXP cleanupG(double *a, double *r, int *e, SEXP beta, SEXP iter, SEXP df, SEXP loss);

SEXP cleanupG_ssr(double *a, double *r, int *e, SEXP beta, SEXP iter, SEXP df, SEXP loss, SEXP rejections) {
  Free(a);
  Free(r);
  Free(e);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, iter);
  SET_VECTOR_ELT(res, 2, df);
  SET_VECTOR_ELT(res, 3, loss);
  SET_VECTOR_ELT(res, 4, rejections);
  UNPROTECT(6);
  return(res);
}

// sequential strong rule
void ssr_glasso(int *screen, double *X, double *r, int *K1, double *lam, double lam_max, int n, int l, int p, int J) {
  
  int K = 0;
  double z_norm, cutoff;

  for (int g = 0; g < J; g++) {
    K = K1[g+1] - K1[g];
    double *z = Calloc(K, double);

    for (int j = K1[g]; j < K1[g+1]; j++) {
      z[j-K1[g]] = crossprod(X, r, n, j) / n;
    }
    z_norm = norm(z, K);
    if (l != 0) {
      cutoff = sqrt(K) * (2 * lam[l] - lam[l-1]);
    } else {
      cutoff = sqrt(K) * (2 * lam[l] - lam_max);
    }

    if (l < 5) {
      Rprintf("\t\t g = %d; K1[g] = %d; z_norm = %f; cutoff = %f, sqrt(%d) = %f\n", g, K1[g], z_norm, cutoff, K, sqrt(K));
    }
    
    if (z_norm > cutoff) {
      screen[g] = 1; // not reject
    } else {
      screen[g] = 0; // reject
    }
    Free(z);
  }
}

// Scan for violations in strong set
int check_strong_set(int *screen, int *e, double *X, double *r, int *K1, double lam, int n, int p, int J) {
  int violations = 0;
  
  int K = 0;
  double z_norm;
  
  for (int g = 0; g < J; g++) {
    if (e[g] == 0 && screen[g] == 1) {
      K = K1[g+1] - K1[g];
      double *z = Calloc(K, double);
      
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      z_norm = norm(z, K);
      if (z_norm > lam * sqrt(K)) {
        e[g] = 1;
        violations++;
      }
      Free(z);
    }
  }
  return violations;
}

// Scan for violations in rest set
int check_rest_set(int *screen, int *e, double *X, double *r, int *K1, double lam, int n, int p, int J) {
  int violations = 0;
  
  int K = 0;
  double z_norm;
  
  for (int g = 0; g < J; g++) {
    if (screen[g] == 0) {
      K = K1[g+1] - K1[g];
      double *z = Calloc(K, double);
      
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      z_norm = norm(z, K);
      if (z_norm > lam * sqrt(K)) {
        e[g] = screen[g] = 1;
        violations++;
      }
      Free(z);
    }
  }
  return violations;
}

SEXP gdfit_gaussian_ssr(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, SEXP lambda, SEXP lam_max_, SEXP alpha_, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, SEXP dfmax_, SEXP gmax_, SEXP user_) {

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
  double lam_max = REAL(lam_max_)[0];
  double alpha = REAL(alpha_)[0];
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double gamma = REAL(gamma_)[0];
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta, iter, df, loss, rejections;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(df = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(df)[i] = 0;
  PROTECT(loss = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(loss)[i] = 0;
  double *b = REAL(beta);
  PROTECT(rejections = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(rejections)[i] = 0;

  // Intermediate quantities
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(J, int);
  for (int g=0; g<J; g++) e[g] = 0;
  int *screen = Calloc(J, int);
  int converged, lstart = 0, ng, nv, violations;
  double shift, l1, l2;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    // lstart = 1;
  }

  Rprintf("J = %d\n", J);
  for (int g = 0; g < J; g++) {
    Rprintf("K1[%d] = %d\n", g, K1[g]);
  }
  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
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
        Free(screen);
        res = cleanupG_ssr(a, r, e, beta, iter, df, loss, rejections);
        return(res);
      }
    }

    // screening
    ssr_glasso(screen, X, r, K1, lam, lam_max, n, l, p, J);
    INTEGER(rejections)[l] = J - sum_rejections(screen, J);
    if (l == 0) Rprintf("lamb_max = %f; \t lam[%d] = %f\n", lam_max, l, lam[l]);
    if (l < 5) {
      Rprintf("screenig: l = %d; lam[l] = %f; rejections[l] = %d\n", l, lam[l], INTEGER(rejections)[l]);
      for (int j=0; j < J; j++) {
        Rprintf("\tscreen[%d] = %d; \te[%d] = %d\n", j, screen[j], j, e[j]);
      }
    }
    
    while (INTEGER(iter)[l] < max_iter) {
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
            if (e[g] && screen[g]) gd_gaussian(b, X, r, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
          }
          
          // Check convergence
          if (checkConvergence(b, a, eps, l, p)) {
            converged  = 1;
            REAL(loss)[l] = gLoss(r,n);
            break;
          }
          for (int j=0; j<p; j++) a[j] = b[l*p+j];
        }
       
        // Scan for violations in strong set
        violations = check_strong_set(screen, e, X, r, K1, lam[l], n, p, J);
        //violations = 0;
        //for (int g=0; g<J; g++) {
        //  if (e[g]==0 && screen[g]) {
        //    l1 = lam[l] * m[g] * alpha;
        //    l2 = lam[l] * m[g] * (1-alpha);
        //    gd_gaussian(b, X, r, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
        //    if (b[l*p+K1[g]] != 0) {
        //      e[g] = 1;
        //      violations++;
        //    }
        //  }
        //}
        if (violations == 0) break;
      }
      
      // Scan for violations in rest set
      violations = check_rest_set(screen, e, X, r, K1, lam[l], n, p, J);
      // Scan for violations
      //violations = 0;
      //for (int g=0; g<J; g++) {
      //  if (screen[g] == 0) {
      //    l1 = lam[l] * m[g] * alpha;
      //    l2 = lam[l] * m[g] * (1-alpha);
      //    gd_gaussian(b, X, r, g, K1, n, l, p, penalty, l1, l2, gamma, df, a);
      //    if (b[l*p+K1[g]] != 0) {
      //      e[g] = screen[g] = 1;
      //      violations++;
      //    }
      //  }
      //}
      if (violations == 0) {
        REAL(loss)[l] = gLoss(r, n);
        break;
      }
    }
  }
  Free(screen);
  res = cleanupG_ssr(a, r, e, beta, iter, df, loss, rejections);
  return(res);
}
