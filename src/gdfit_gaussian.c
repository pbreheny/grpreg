#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double crossprod(double *x, double *y, int n, int j);
double norm(double *x, int p);
double S(double, double);
double F(double, double, double, double);
double Fs(double, double, double, double);
double gLoss(double *r, int n);

// Sum rejections
int sum_rejections(int *x, int n) {
  double val = 0;
  for (int i=0; i<n; i++) val += x[i];
  return(val);
}

// Group descent update
void gd_gaussian(double *b, double *x, double *r, int g, int *K1, int *K,
                       int n, int l, int p, const char *penalty, double lam1, double lam2,
                       double gamma, SEXP df, double *a, double *maxChange) {
  // Calculate z
  double *z = Calloc(K[g], double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + a[j];
  double z_norm = norm(z, K[g]);
  
  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(z_norm, lam1) / (1+lam2);
  if (strcmp(penalty, "grMCP")==0) len = F(z_norm, lam1, lam2, gamma);
  if (strcmp(penalty, "grSCAD")==0) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 || a[K1[g]] != 0) {
    // If necessary, update beta and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      b[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = b[l*p+j]-a[j];
      if (fabs(shift) > maxChange[0]) maxChange[0] = fabs(shift);
      for (int i=0; i<n; i++) r[i] -= x[n*j+i] * shift;
    }
  }
  // Update df
  if (len > 0) REAL(df)[l] += K[g] * len / z_norm;
  Free(z);
}

// Scan for violations in rest set
int check_rest_set(int *e2, int *e, double *xTr, double *X, double *r, int *K1, int *K, double lam, int n, int J) {
  int violations = 0;
  double TOLERANCE = 1e-8;
  for (int g = 0; g < J; g++) {
    if (e2[g] == 0) {
      double *z = Calloc(K[g], double);
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      xTr[g] = norm(z, K[g]);
      if (xTr[g] + TOLERANCE > lam * sqrt(K[g])) {
        e[g] = e2[g] = 1;
        violations++;
      }
      Free(z);
    }
  }
  return violations;
}

// Scan for violations in strong set
int check_strong_set(int *e2, int *e, double *xTr, double *X, double *r, int *K1, int *K, double lam, int n, int J) {
  int violations = 0;
  for (int g = 0; g < J; g++) {
    if (e[g] == 0 && e2[g] == 1) {
      double *z = Calloc(K[g], double);
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      xTr[g] = norm(z, K[g]);
      if (xTr[g] > lam * sqrt(K[g])) {
        e[g] = 1;
        violations++;
      }
      Free(z);
    }
  }
  return violations;
}

// sequential strong rule
void ssr_glasso(int *e2, double *xTr, int *K1, int *K, double *lam, double lam_max, int l, int J) {
  double cutoff;
  double TOLERANCE = 1e-8;
  for (int g = 0; g < J; g++) {
    if (l != 0) {
      cutoff = sqrt(K[g]) * (2 * lam[l] - lam[l-1]);
    } else {
      if (lam_max > 0) cutoff = sqrt(K[g]) * (2 * lam[l] - lam_max);
      else cutoff = 0;
    }
    if (xTr[g] + TOLERANCE > cutoff) {
      e2[g] = 1; // not reject, in strong set
    } else {
      e2[g] = 0; // reject
    }
  }
}

// BEDPP initialization
void bedpp_init(double *yTxxTv1, double *xTv1_sq, double *xTy_sq, double *xTr,
                double *X, double *y, int *K1, int *K, int *g_star_ptr, 
                int *K_star_ptr, int K1_len, int n, int J) {

  // compute Xj^T * y
  double tmp = 0;
  double *XTy = Calloc(K1_len, double);
  for (int g=0; g<J; g++) {
    for (int j = K1[g]; j < K1[g+1]; j++) {
      XTy[j-K1[0]] = crossprod(X, y, n, j); // K1 contains consecutive indices of penalized groups
      xTy_sq[g] += pow(XTy[j-K1[0]], 2);
    }
    xTr[g] = sqrt(xTy_sq[g]) / n;
    if (xTr[g] / sqrt(K[g]) > tmp) {
      tmp = xTr[g] / sqrt(K[g]);
      *g_star_ptr = g;
      *K_star_ptr = K[g];
    }
  }

  // compute v1_bar at lam_max: = X* X*^T y
  double *v1_lam_max_tmp = Calloc(*K_star_ptr, double);
  double *v1_lam_max = Calloc(n, double); // tmp quantity for BEDPP: v1_bar at lam_max: = X_star * X_star^T * y
  for (int j = K1[*g_star_ptr]; j < K1[*g_star_ptr+1]; j++) {
    v1_lam_max_tmp[j-K1[*g_star_ptr]] = crossprod(X, y, n, j);
  }
  for (int i = 0; i < n; i++) {
    for (int j = K1[*g_star_ptr]; j < K1[*g_star_ptr+1]; j++) {
      v1_lam_max[i] += X[i+n*j] * v1_lam_max_tmp[j-K1[*g_star_ptr]];
    }
  }
  
  // compute Xj^T * v1, then compute its norm_sq, then compute  y^T * Xj * Xj^T * v1
  double xTv1_tmp;
  for (int g = 0; g < J; g++) {
    for (int j = K1[g]; j < K1[g+1]; j++) {
      xTv1_tmp = crossprod(X, v1_lam_max, n, j);
      xTv1_sq[g] += pow(xTv1_tmp, 2);
      yTxxTv1[g] += XTy[j-K1[0]] * xTv1_tmp; 
    }
  }
  Free(XTy);
  Free(v1_lam_max_tmp);
  Free(v1_lam_max);
}

// basic EDPP screening
void bedpp_glasso(int *e3, double *yTxxTv1, double *xTv1_sq, double *xTy_sq, 
                  double y_norm_sq, int *K, double lam, double lam_max, 
                  int K_star, int n, int J) {
  double TOLERANCE = 1e-8;
  double LHS, RHS, LHS_temp;
  for (int g = 0; g < J; g++) {
    LHS_temp = pow(lam + lam_max, 2) * xTy_sq[g] - 2 * (lam_max * lam_max - lam * lam) * yTxxTv1[g] / n +
      pow((lam_max - lam) / n, 2) * xTv1_sq[g];
    if (LHS_temp < 0) {
      LHS = 0.0;
    } else {
      LHS = sqrt(LHS_temp);
    }
    RHS = 2 * n * lam * lam_max * sqrt(K[g]) - (lam_max - lam) * sqrt(n * y_norm_sq - pow(n * lam_max, 2) * K_star);
    
    if (LHS + TOLERANCE > RHS) {
      e3[g] = 1; // not reject, thus in BEDPP set
    } else {
      e3[g] = 0; // reject
    }
  }
}

// hybrid sequential safe-strong rule: SSR-BEDPP
void ssr_bedpp_glasso(int *e2, int *e3, double *xTr, int *K1, int *K, double *lam, double lam_max, int l, int J) {
  double cutoff;
  for (int g = 0; g < J; g++) {
    if (e3[g] == 1) { // only check SSR in BEDPP set
      if (l != 0) {
        cutoff = sqrt(K[g]) * (2 * lam[l] - lam[l-1]);
      } else {
        cutoff = sqrt(K[g]) * (2 * lam[l] - lam_max);
      }
      if (xTr[g] > cutoff) {
        e2[g] = 1; // not reject
      } else {
        e2[g] = 0; // reject
      }
    } else {
      e2[g] = 0; // reject
    }
  }
}

// update xTr[g] for groups which are rejected at previous lambda but accepted at current one.
void update_xTr(int *e3, int *e3_old, double *xTr, double *X, double *r, int *K1, int *K, int n, int J) {
  for (int g = 0; g < J; g++) {
    if (e3[g] == 1 && e3_old[g] == 0) {
      double *z = Calloc(K[g], double);
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      xTr[g] = norm(z, K[g]);
    }
  }
}

// check rest set with SSR-BEDPP screening
int check_rest_set_ssr_bedpp(int *e2, int *e, int *e3, double *xTr, double *X, 
                             double *r, int *K1, int *K, double lam, int n, int J) {
  int violations = 0;
  for (int g = 0; g < J; g++) {
    if (e3[g] == 1 && e2[g] == 0) { // check groups not rejected by BEDPP but by SSR
      double *z = Calloc(K[g], double);
      for (int j = K1[g]; j < K1[g+1]; j++) {
        z[j-K1[g]] = crossprod(X, r, n, j) / n;
      }
      xTr[g] = norm(z, K[g]);
      if (xTr[g] > lam * sqrt(K[g])) {
        e[g] = e2[g] = 1;
        violations++;
      }
      Free(z);
    }
  }
  return violations;
}


SEXP gdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP K1_, SEXP K0_, 
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
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(group_multiplier);
  int dfmax = INTEGER(dfmax_)[0];
  int gmax = INTEGER(gmax_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta, iter, df, loss, rejections, safe_rejections;
  PROTECT(res = allocVector(VECSXP, 6));
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
  PROTECT(safe_rejections = allocVector(INTSXP, L)); // # of rejections by BEDPP
  for (int i=0; i<L; i++) INTEGER(safe_rejections)[i] = 0;

  // Intermediate quantities
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j] = 0;
  int *e = Calloc(J, int); // ever-active set
  for (int g=0; g<J; g++) e[g] = 0;
  int converged, lstart = 0, ng, nv, violations;
  double shift, l1, l2, maxChange;

  // variables for screening
  int *e2 = Calloc(J, int); // strong set
  int *e3 = Calloc(J, int); // BEDPP set
  int *e3_old = Calloc(J, int); // previous BEDPP set
  for (int g=0; g<J; g++) e3_old[g] = 1; // initialize e3_old = 1, not reject.
  int *K = Calloc(J, int); // group size
  int K1_len = 0; // # of variables in K1
  for (int g = 0; g < J; g++) {
    K[g] = K1[g+1] - K1[g];
    K1_len += K[g];
  }
  double *xTr = Calloc(J, double);
  double *xTy_sq = Calloc(J, double); // tmp quantity for BEDPP: square norm of X^T*y
  double *yTxxTv1 = Calloc(J, double); // tmp quantity for BEDPP: y^T*X*X^T v1
  double *xTv1_sq = Calloc(J, double); // tmp quantity for BEDPP: square norm of X^T*v1
  int g_star, K_star; // group index and size corresponding to lambda_max
  int *g_star_ptr = &g_star;
  int *K_star_ptr = &K_star;
  double y_norm_sq = pow(norm(y, n), 2);
  int bedpp_flag;
  if ((strcmp(penalty, "grLasso")==0) & !user) {
    bedpp_flag = 1;
    bedpp_init(yTxxTv1, xTv1_sq, xTy_sq, xTr, X, y, K1, K, g_star_ptr, K_star_ptr, K1_len, n, J);
  }
  else bedpp_flag = 0;
  
  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  double rss = gLoss(r,n);
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = rss;
    INTEGER(rejections)[0] = J;
    INTEGER(safe_rejections)[0] = J;
    lstart = 1;
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
    if (bedpp_flag) { //SSR-BEDPP screening
      // BEDPP screening
      bedpp_glasso(e3, yTxxTv1, xTv1_sq, xTy_sq, y_norm_sq, K, lam[l], lam_max, K_star, n, J);
      INTEGER(safe_rejections)[l] = J - sum_rejections(e3, J);
      // update xTr[g] for groups which are rejected at previous lambda but accepted at current one.
      update_xTr(e3, e3_old, xTr, X, r, K1, K, n, J);
      for (int g=0; g<J; g++) e3_old[g] = e3[g]; // reset e3_old to be new e3;
      // SSR screening
      ssr_bedpp_glasso(e2, e3, xTr, K1, K, lam, lam_max, l, J);
    } else { // only SSR screening
      INTEGER(safe_rejections)[l] = 0;
      ssr_glasso(e2, xTr, K1, K, lam, lam_max, l, J);
    }
    INTEGER(rejections)[l] = J - sum_rejections(e2, J);
    if (INTEGER(safe_rejections)[l] <= 0) bedpp_flag = 0; // BEDPP not effective, turn it off for next lambda.
    
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {
          converged = 0;
          INTEGER(iter)[l]++;
          tot_iter++;
          REAL(df)[l] = 0;
          maxChange = 0;
          
          // Update unpenalized covariates
          for (int j=0; j<K0; j++) {
            shift = crossprod(X, r, n, j)/n;
            if (fabs(shift) > maxChange) maxChange = fabs(shift);
            b[l*p+j] = shift + a[j];
            for (int i=0; i<n; i++) r[i] -= shift * X[n*j+i];
            REAL(df)[l] += 1;
          }
          
          // Update penalized groups
          for (int g=0; g<J; g++) {
            l1 = lam[l] * m[g] * alpha;
            l2 = lam[l] * m[g] * (1-alpha);
            if (e[g]) {
              gd_gaussian(b, X, r, g, K1, K, n, l, p, penalty, l1, l2, gamma, df, a, &maxChange);
            }
          }
          
          // Check convergence
	  for (int j=0; j<p; j++) a[j] = b[l*p+j];
	  if (maxChange < eps*sdy) break;
        }
        // Scan for violations in strong set
        violations = check_strong_set(e2, e, xTr, X, r, K1, K, lam[l], n, J);
        if (violations == 0) break;
      }

      // Scan for violations in rest set
      if (bedpp_flag) {
        violations = check_rest_set_ssr_bedpp(e2, e, e3, xTr, X, r, K1, K, lam[l], n, J);
      } else {
        violations = check_rest_set(e2, e, xTr, X, r, K1, K, lam[l], n, J);
      }
      if (violations == 0) {
        REAL(loss)[l] = gLoss(r, n);
        break;
      }
    }
  }
  Free(a);
  Free(r);
  Free(e);
  Free(e2);
  Free(e3);
  Free(e3_old);
  Free(K);
  Free(yTxxTv1);
  Free(xTv1_sq);
  Free(xTy_sq);
  Free(xTr);
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, iter);
  SET_VECTOR_ELT(res, 2, df);
  SET_VECTOR_ELT(res, 3, loss);
  SET_VECTOR_ELT(res, 4, rejections);
  SET_VECTOR_ELT(res, 5, safe_rejections);
  UNPROTECT(7);
  return(res);
}
