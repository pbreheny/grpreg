#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double norm(double *x, int p);

SEXP maxprod(SEXP X_, SEXP y_, SEXP K_, SEXP m_) {

  // Declarations
  int n = nrows(X_);
  int J = length(K_)-1;
  SEXP zmax;
  PROTECT(zmax = allocVector(REALSXP, 1));
  REAL(zmax)[0] = 0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *K = INTEGER(K_);

  for (int g=0; g<J; g++) {
    for (int j=K[g]; j<K[g+1]; j++) {
      double z = fabs(crossprod(X, y, n, j) / m[g]);
      if (z > REAL(zmax)[0]) REAL(zmax)[0] = z;
    }
  }

  // Return
  UNPROTECT(1);
  return(zmax);
}

SEXP maxgrad(SEXP X_, SEXP y_, SEXP K_, SEXP m_) {

  // Declarations
  int n = nrows(X_);
  int J = length(K_)-1;
  SEXP zmax;
  PROTECT(zmax = allocVector(REALSXP, 1));
  REAL(zmax)[0] = 0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *K = INTEGER(K_);

  for (int g=0; g<J; g++) {
    int Kg = K[g+1]-K[g];
    double *Z = Calloc(Kg, double);
    for (int j=K[g]; j<K[g+1]; j++) {
      Z[j-K[g]] = crossprod(X, y, n, j);
    }
    double z = norm(Z, Kg) / m[g];
    if (z > REAL(zmax)[0]) REAL(zmax)[0] = z;
    Free(Z);
  }

  // Return
  UNPROTECT(1);
  return(zmax);
}
