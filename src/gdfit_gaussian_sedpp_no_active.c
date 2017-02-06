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
double gLoss(double *r, int n);
int sum_rejections(int *x, int n);
// Group descent update
void gd_gaussian_ssr(double *b, double *x, double *r, int g, int *K1, int *K, int n, int l, int p, const char *penalty, double lam1, double lam2, double gamma, SEXP df, double *a);
// pv2 = V2 - <v1, v2> / ||v1||^2_2 * V1
void update_pv2(double *pv2, double *v1, double *v2, int n);
// sequential EDPP rule
void sedpp_glasso(int *e2, double *X, double *r, double *y, double *v1_bar_lam_max, int *K1, double *lam, double lam_max, int n, int p, int l, int J);

SEXP cleanupG_sedpp_no_active(double *a, double *r, int *e2, int *K, double *xTr, 
                    double *v1_bar_lam_max, double *v1_bar_lam_max_tmp,
                    SEXP beta, SEXP iter, SEXP df, SEXP loss, SEXP rejections) {
  Free(a);
  Free(r);
  Free(e2);
  Free(K);
  Free(xTr);
  Free(v1_bar_lam_max);
  Free(v1_bar_lam_max_tmp);
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

SEXP gdfit_gaussian_sedpp_no_active(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, 
                        SEXP lambda, SEXP lam_max_, SEXP alpha_, SEXP eps_, 
                        SEXP max_iter_, SEXP gamma_, SEXP group_multiplier, 
                        SEXP dfmax_, SEXP gmax_, SEXP user_) {
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
  int *K = Calloc(J, int); // group size
  int converged, lstart = 0, ng, nv, violations;
  double shift, l1, l2;

  // variables for screening
  int g_star; // group index corresponding to lambda_max
  int *e2 = Calloc(J, int); // EDPP set
  double tmp = 0;
  double *xTr = Calloc(J, double);
  for (int g=0; g<J; g++) {
    K[g] = K1[g+1] - K1[g];
    double *z = Calloc(K[g], double);
    for (int j = K1[g]; j < K1[g+1]; j++) {
      z[j-K1[g]] = crossprod(X, r, n, j) / n;
    }
    xTr[g] = norm(z, K[g]);
    if (xTr[g] / sqrt(K[g]) > tmp) {
      tmp = xTr[g] / sqrt(K[g]) ;
      g_star = g;
    }
    Free(z);
  }

  // compute v1_bar at lam_max: = X* X*^T y
  int K_star = K1[g_star+1] - K1[g_star];
  double *v1_bar_lam_max = Calloc(n, double);
  for (int i=0; i<n; i++) v1_bar_lam_max[i] = 0;
  double *v1_bar_lam_max_tmp = Calloc(K_star, double);
  for (int j = K1[g_star]; j < K1[g_star+1]; j++) {
    v1_bar_lam_max_tmp[j-K1[g_star]] = crossprod(X, y, n, j);
  }
  for (int i = 0; i < n; i++) {
    for (int j = K1[g_star]; j < K1[g_star+1]; j++) {
      v1_bar_lam_max[i] += X[i+n*j] * v1_bar_lam_max_tmp[j-K1[g_star]];
    }
  }

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    INTEGER(rejections)[0] = J;
    lstart = 1;
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
        res = cleanupG_sedpp_no_active(a, r, e2, K, xTr, v1_bar_lam_max, 
                                       v1_bar_lam_max_tmp, beta, iter, df, loss, rejections);
        return(res);
      }
    }

    // screening
    sedpp_glasso(e2, X, r, y, v1_bar_lam_max, K1, lam, lam_max, n, p, l, J);
    INTEGER(rejections)[l] = J - sum_rejections(e2, J);
    
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
        if (e2[g]) {
          gd_gaussian_ssr(b, X, r, g, K1, K, n, l, p, penalty, l1, l2, gamma, df, a);
        }
      }
      
      // Check convergence
      if (checkConvergence(b, a, eps, l, p)) {
        converged  = 1;
        REAL(loss)[l] = gLoss(r,n);
        break;
      }
      for (int j=0; j<p; j++) a[j] = b[l*p+j];
    }
  }
  
  res = cleanupG_sedpp_no_active(a, r, e2, K, xTr, v1_bar_lam_max, v1_bar_lam_max_tmp,
                       beta, iter, df, loss, rejections);
  return(res);
}
