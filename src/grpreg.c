#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

// Allocates memory for a vector of size n
static double *vector(int n) {
  double *v;
  v = Calloc(n, double);
  return v;
}

// Check for convergence of beta[l]
static int checkConvergence(double *beta, double *beta_old, double eps, int l, int J)
{
  int j;
  int converged = 1;
  for (j=0; j < J; j++) {
    if (beta[l*J+j]!=0 & beta_old[j]!=0) {
      if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps) {
	converged = 0;
	break;
      }
    } else if (beta[l*J+j]==0 & beta_old[j]!=0) {
      converged = 0;
      break;
    } else if (beta[l*J+j]!=0 & beta_old[j]==0) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

// Soft-thresholding operator
static double S(double z, double l)
{
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

// Firm-thresholding operator
static double F(double z, double l1, double l2, double gamma)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(1+l2-1/gamma));
  else return(z/(1+l2));
}

// SCAD-modified firm-thresholding operator
static double Fs(double z, double l1, double l2, double gamma)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(1+l2));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(1-1/(gamma-1)+l2));
  else return(z/(1+l2));
}

// Firm-thresholding operator w/ adaptive rescaling
static double F_ars(double z, double l1, double l2, double gamma, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}

// SCAD-modified firm-thresholding operator
static double Fs_ars(double z, double l1, double l2, double gamma, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}

// MCP penalty
static double MCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta <= a*l) return(l*theta - pow(theta,2)/(2*a));
  else return(a*pow(l,2)/2);
}

// MCP penalization rate
static double dMCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta < a*l) return(l-theta/a);
  else return(0);
}

// Group descent update
static void gd_gaussian(double *beta, double *x, double *r, int K0, int Kj, int n, int l, int p, char *penalty, double lam1, double lam2, double gamma, double *df, double *beta_old)
{
  // Calculate z
  int K = Kj - K0;
  double *z;
  z = vector(K);
  for (int j=K0; j<Kj; j++) {
    for (int i=0; i<n; i++) z[j-K0] = z[j-K0] + x[n*j+i]*r[i];
    z[j-K0] = z[j-K0]/n + beta_old[j];
  }
  double z_norm = 0;
  for (int j=0; j<K; j++) z_norm = z_norm + pow(z[j],2);
  z_norm = sqrt(z_norm);

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(z_norm, lam1) / (1+lam2);
  if (strcmp(penalty, "grMCP")==0) len = F(z_norm, lam1, lam2, gamma);
  if (strcmp(penalty, "grSCAD")==0) len = Fs(z_norm, lam1, lam2, gamma);
  for (int j=K0; j<Kj; j++) beta[l*p+j] = len * z[j-K0] / z_norm;
  //Rprintf("znorm: %f   len: %f\n", z_norm, len);

  // Update r
  for (int j=K0; j<Kj; j++) {
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
    }
  }

  // Update df
  df[0] = df[0] + K * len / z_norm;
  Free(z);
}

// Group descent update -- binomial
static void gd_binomial(double *beta, double *x, double *r, int K0, int Kj, int n, int l, int p, char *penalty, double lam1, double lam2, double gamma, double *df, double *beta_old, double *w)
{
  // Calculate z, v
  int K = Kj - K0;
  double *z, *v;
  z = vector(K);
  v = vector(K);
  for (int j=K0; j<Kj; j++) {
    for (int i=0; i<n; i++) {
      z[j-K0] = z[j-K0] + x[n*j+i] * w[i] * r[i];
      v[j-K0] = v[j-K0] + x[n*j+i] * w[i] * x[n*j+i];
    }
    v[j-K0] = v[j-K0]/n;
    z[j-K0] = z[j-K0]/n + v[j-K0]*beta_old[j];
  }
  double z_norm = 0;
  for (int j=0; j<K; j++) z_norm = z_norm + pow(z[j],2);
  z_norm = sqrt(z_norm);
  double vbar = 0;
  for (int j=0; j<K; j++) vbar = vbar + v[j];
  vbar = vbar / n;

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(z_norm, lam1);
  if (strcmp(penalty, "grMCP")==0) len = F(z_norm, lam1, lam2, gamma);
  if (strcmp(penalty, "grSCAD")==0) len = Fs(z_norm, lam1, lam2, gamma);
  for (int j=K0; j<Kj; j++) beta[l*p+j] = (len * z[j-K0]) / (z_norm * v[j-K0] * (1+lam2));
  //Rprintf("znorm: %f   len: %f\n", z_norm, len);

  // Update r
  for (int j=K0; j<Kj; j++) {
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
    }
  }

  // Update df
  df[0] = df[0] + K * len / z_norm;
  Free(z);
  Free(v);
}

// Groupwise local coordinate descent updates
static void gLCD_gaussian(double *beta, char *penalty, double *x, double *r, int K0, int Kj, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *beta_old)
{
  // Make initial local approximation
  int K = Kj - K0;
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "geLasso")==0) for (int j=K0; j<Kj; j++) sG = sG + fabs(beta[l+p+j]);
  if (strcmp(penalty, "geMCP")==0) for (int j=K0; j<Kj; j++) sG = sG + MCP(beta[l+p+j], lam1, gamma);
  if (strcmp(penalty, "gMCP")==0) for (int j=K0; j<Kj; j++) sG = sG + MCP(beta[l+p+j], lam1, gamma);

  // Coordinate descent
  for (int j=K0; j<Kj; j++) {

    // Calculate LS sol'n
    double z=0;
    for (int i=0; i<n; i++) z = z + x[n*j+i]*r[i];
    z = z/n + beta_old[j];

    // Calculate ljk
    double ljk=0;
    if (lam1 != 0) {
      if (strcmp(penalty, "gMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(beta[l*p+j], lam1, gamma);
      if (strcmp(penalty, "geLasso")==0) ljk = lam1*exp(-tau/lam1*sG);
      if (strcmp(penalty, "geMCP")==0) ljk = exp(-tau/pow(lam1,2)*sG)*dMCP(beta[j],lam1,gamma);
    }

    // Update beta
    beta[l*p+j] = S(z, ljk) / (1+lam2);
    //Rprintf("z: %f  ljk:%f  b:%f  oldbeta: %f\n", z, ljk, beta[l*p+j], beta_old[j]);

    // Update r
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
      if (strcmp(penalty, "geLasso")==0) sG = sG + fabs(beta[l+p+j]) - fabs(beta_old[j]);
      if (strcmp(penalty, "geMCP")==0) sG = sG + MCP(beta[l+p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
      if (strcmp(penalty, "gMCP")==0) sG = sG + MCP(beta[l+p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
    }

    // Update df
    df[0] = df[0] + fabs(beta[l*j+p])/fabs(z);
    }
}

// Groupwise local coordinate descent updates -- binomial
static void gLCD_binomial(double *beta, char *penalty, double *x, double *r, int K0, int Kj, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *beta_old, double *w)
{
  // Calculate v
  int K = Kj - K0;
  double *v;
  v = vector(K);
  for (int j=K0; j<Kj; j++) {
    for (int i=0; i<n; i++) v[j-K0] = v[j-K0] + x[j*n+i] * w[i] * x[j*n+i];
    v[j-K0] = v[j-K0]/n;
  }

  // Make initial local approximation
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "geLasso")==0) for (int j=K0; j<Kj; j++) sG = sG + fabs(v[j-K0] * beta[l+p+j]);
  if (strcmp(penalty, "geMCP")==0) for (int j=K0; j<Kj; j++) sG = sG + MCP(v[j-K0] * beta[l+p+j], lam1, gamma);
  if (strcmp(penalty, "gMCP")==0) for (int j=K0; j<Kj; j++) sG = sG + MCP(v[j-K0] * beta[l+p+j], lam1, gamma);

  // Coordinate descent
  for (int j=K0; j<Kj; j++) {

    // Calculate z
    double z = 0;
    for (int i=0; i<n; i++) z = z + x[j*n+i] * w[i] * r[i];
    z = z/n + v[j-K0]*beta_old[j];

    // Calculate ljk
    double ljk=0;
    if (lam1 != 0) {
      if (strcmp(penalty, "gMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(beta[l*p+j], lam1, gamma);
      if (strcmp(penalty, "geLasso")==0) ljk = lam1*exp(-tau/lam1*sG);
      if (strcmp(penalty, "geMCP")==0) ljk = exp(-tau/pow(lam1,2)*sG)*dMCP(beta[j],lam1,gamma);
    }

    // Update beta
    beta[l*p+j] = S(z, ljk) / (v[j-K0]*(1+lam2));
    //Rprintf("z: %f  ljk:%f  b:%f  oldbeta: %f\n", z, ljk, beta[l*p+j], beta_old[j]);

    // Update r
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
      if (strcmp(penalty, "geLasso")==0) sG = sG + fabs(v[j-K0]*beta[l+p+j]) - fabs(v[j-K0]*beta_old[j]);
      if (strcmp(penalty, "geMCP")==0) sG = sG + MCP(v[j-K0]*beta[l+p+j], lam1, gamma) - MCP(v[j-K0]*beta_old[j], lam1, gamma);
      if (strcmp(penalty, "gMCP")==0) sG = sG + MCP(v[j-K0]*beta[l+p+j], lam1, gamma) - MCP(v[j-K0]*beta_old[j], lam1, gamma);
    }

    // Update df
    df[0] = df[0] + fabs(beta[l*j+p])/fabs(z);
    }
  Free(v);
}

static void gpPathFit_binomial(double *beta0, double *beta, int *iter, double *df, double *x, double *y, int *group, int *n_, int *p_, char **penalty_, int *J_, int *K, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *tau_, double *group_multiplier, int *dfmax_, int *warn_)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; int warn = warn_[0];
  double beta0_old;
  double *r, *beta_old, *w;
  r = vector(n);
  w = vector(n);
  beta_old = vector(p);

  // Initialization
  double ybar=0;
  for (int i=0; i<n; i++) ybar = ybar + y[i];
  ybar = ybar/n;
  double nullDev = 0;
  for (int i=0;i<n;i++) nullDev = nullDev - y[i]*log(ybar) - (1-y[i])*log(1-ybar);

  // Path
  double xwr, xwx, eta, pi, z, v, Dev;
  for (int l=0; l<L; l++) {
    if (l != 0) {
      beta0_old = beta0[l-1];
      for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    }
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l] = iter[l] + 1;

      // Check dfmax
      int active = 0;
      for (int j=0; j<p; j++) if (beta[l*p+j]!=0) active++;
      if (active > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	return;
      }

      // Approximate L
      Dev = 0;
      for (int i=0; i<n; i++) {
	eta = beta0_old;
	for (int j=0; j<p; j++) eta = eta + x[j*n+i] * beta_old[j];
	pi = exp(eta) / (1+exp(eta));
	if (pi > .9999) {
	  pi = 1;
	  w[i] = .0001;
	} else if (pi < .0001) {
	  pi = 0;
	  w[i] = .0001;
	} else w[i] = pi*(1-pi);

	r[i] = (y[i] - pi) / w[i];
	Dev = Dev - y[i]*log(pi) - (1-y[i])*log(1-pi);
      }

      // Check for saturation
      if (Dev/nullDev < .01) {
	if (warn) warning("Model saturated; exiting...");
	for (int ll=l;ll<L;ll++) {
	  beta0[ll] = R_NaReal;
	  for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(w);
	Free(r);
	return;
      }
      df[l] = 0;

      // Update intercept
      xwr = xwx = 0;
      for (int i=0;i<n;i++) {
	xwr = xwr + w[i]*r[i];
	xwx = xwx + w[i];
      }
      beta0[l] = xwr/xwx + beta0_old;
      for (int i=0;i<n;i++) r[i] = r[i] - (beta0[l] - beta0_old);
  
      // Update unpenalized covariates
      int j;
      for (j=0;; j++) {
	if (group[j]!=0) break;
	xwr=0;
	xwx=0;
	for (int i=0; i<n; i++) {
	  xwr = xwr + x[j*n+i]*w[i]*r[i];
	  xwx = xwx + x[j*n+i]*w[i]*x[j*n+i];
	}
	z = xwr/n;
	v = xwx/n;
	beta[l*p+j] = z/v + beta_old[j];
        for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[j*n+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	int K0 = j;
	int Kj = j + K[g];
	if (strncmp(penalty, "gr", 2)==0) gd_binomial(beta, x, r, K0, Kj, n, l, p, penalty, lam1[l]*group_multiplier[g], lam2[l], gamma, &df[l], beta_old, w);
	else gLCD_binomial(beta, penalty, x, r, K0, Kj, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, &df[l], beta_old, w);
	j = Kj;
      }

      // Check convergence
      if (checkConvergence(beta, beta_old, eps, l, p)) {
	converged  = 1;
	break;
      }
      beta0_old = beta0[l];
      for (int j=0; j<p; j++) beta_old[j] = beta[l*p+j];
    }
  }
  Free(beta_old);
  Free(w);
  Free(r);
}

static void gpPathFit_gaussian(double *beta, int *iter, double *df, double *x, double *y, int *group, int *n_, int *p_, char **penalty_, int *J_, int *K, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *tau_, int *dfmax_, int *group_multiplier)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; 
  double *r, *beta_old;
  r = vector(n);
  for (int i=0; i<n; i++) r[i] = y[i];
  beta_old = vector(p);
  //for (int j=0; j < *p; j++) beta_old[j] = 0;

  /* Path */
  for (int l=1; l<L; l++) {
    //Rprintf("lambda: %f %f\n",lam1[l], lam2[l]);
    if (l != 0) for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l] = iter[l] + 1;

      /* Check dfmax */
      int active = 0;
      for (int j=0; j<p; j++) if (beta[l*p+j]!=0) active++;
      if (active > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	return;
      }
      df[l] = 0;

      /* Update unpenalized covariates */
      //Rprintf("b: %f  r1:%f  r2:%f\n", beta[l*p+0],  r[0],  r[1]);
      int j;
      for (j=0;; j++) {
	if (group[j]!=0) break;
	double z=0;
	for (int i=0; i<n; i++) z = z + x[j*n+i]*r[i];
	beta[l*p+j] = z/n + beta_old[j];
        for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[j*n+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	int K0 = j;
	int Kj = j + K[g];
	if (strncmp(penalty, "gr", 2)==0) gd_gaussian(beta, x, r, K0, Kj, n, l, p, penalty, lam1[l]*group_multiplier[g], lam2[l], gamma, &df[l], beta_old);
	else gLCD_gaussian(beta, penalty, x, r, K0, Kj, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, &df[l], beta_old);
	j = Kj;
      }

      // Check for convergence      
      if (checkConvergence(beta, beta_old, eps, l, p)) {
	converged  = 1;
	break;
       }
      for (j=0; j<p; j++) beta_old[j] = beta[l*p+j];
    }
  }
  Free(beta_old);
  Free(r);
}

static const R_CMethodDef cMethods[] = {
  {"gpPathFit_gaussian", (DL_FUNC) &gpPathFit_gaussian, 20},
  {"gpPathFit_binomial", (DL_FUNC) &gpPathFit_binomial, 22},
  NULL
};

void R_init_grpreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
