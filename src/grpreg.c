#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

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

// Cross product of the jth column of x with y
static double crossprod(double *x, double *y, int n, int j)
{
  double val = 0;
  int nn = n*j;
  for (int i=0; i<n; i++) val += x[nn+i] * y[i];
  return(val);
}

// Mean of x
static double mean(double *x, int n)
{
  double val = 0;
  for (int i=0; i<n; i++) val += x[i];
  return(val/n);
}

// Gaussian loss (squared error loss)
static double gLoss(double *r, int n)
{
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Euclidean norm
static double norm(double *x, int p)
{
  double x_norm = 0;
  for (int j=0; j<p; j++) x_norm = x_norm + pow(x[j],2);
  x_norm = sqrt(x_norm);
  return(x_norm);
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
static void gd_gaussian(double *beta, double *x, double *r, int g, int *K1, int *active, int n, int l, int p, char *penalty, double lam1, double lam2, double gamma, double *df, double *beta_old)
{
  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + beta_old[j];
  double z_norm = norm(z,K);

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(z_norm, lam1) / (1+lam2);
  if (strcmp(penalty, "grMCP")==0) len = F(z_norm, lam1, lam2, gamma);
  if (strcmp(penalty, "grSCAD")==0) len = Fs(z_norm, lam1, lam2, gamma);
  if (len != 0 | active[g]) {
    // If necessary, update beta and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      beta[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = beta[l*p+j]-beta_old[j];
      for (int i=0; i<n; i++) r[i] -= x[n*j+i] * shift;
    }
    if (len > 0) active[g] = 1; else active[g] = 0;
  }

  // Update df
  if (len > 0) *df = *df + K * len / z_norm;
  Free(z);
}

// Group descent update -- binomial
static void gd_binomial(double *beta, double *x, double *r, int g, int *K1, int *active, int n, int l, int p, char *penalty, double lam1, double lam2, double gamma, double *df, double *beta_old)
{
  // Calculate z
  int K = K1[g+1] - K1[g];
  double *z = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) z[j-K1[g]] = crossprod(x, r, n, j)/n + beta_old[j];
  double z_norm = norm(z,K);
  double v = 0.25;

  // Update b
  double len;
  if (strcmp(penalty, "grLasso")==0) len = S(v * z_norm, lam1) / (v * (1 + lam2));
  if (strcmp(penalty, "grMCP")==0) len = F(v * z_norm, lam1, lam2, gamma) / v;
  if (strcmp(penalty, "grSCAD")==0) len = Fs(v * z_norm, lam1, lam2, gamma) / v;
  if (len != 0 | active[g]) {
    // If necessary, update beta and r
    for (int j=K1[g]; j<K1[g+1]; j++) {
      beta[l*p+j] = len * z[j-K1[g]] / z_norm;
      double shift = beta[l*p+j]-beta_old[j];
      for (int i=0; i<n; i++) r[i] -= x[n*j+i] * shift;
    }
    if (len > 0) active[g] = 1; else active[g] = 0;
  }

  // Update df
  df[0] = df[0] + K * len / z_norm;
  Free(z);
}

// Groupwise local coordinate descent updates
static void gLCD_gaussian(double *beta, char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *beta_old, double delta)
{
  // Make initial local approximation
  int K = K1[g+1] - K1[g];
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "geLasso")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(beta_old[j]);
  if (strcmp(penalty, "geMCP")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(beta_old[j], lam1, gamma);
  if (strcmp(penalty, "gMCP")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(beta_old[j], lam1, gamma);
  if (strcmp(penalty, "gBridge")==0) {
      for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(beta_old[j]);
      if (sG==0) return;
      if (sG < delta) {
	for (int j=K1[g]; j<K1[g+1]; j++) {
	  beta[l*p+j] = 0;
	  for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
	}
	return;
    }
  }

  // Coordinate descent
  for (int j=K1[g]; j<K1[g+1]; j++) {

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
      if (strcmp(penalty,"gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
    }

    // Update beta
    beta[l*p+j] = S(z, ljk) / (1+lam2);

    // Update r
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
      if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(beta[l*p+j]) - fabs(beta_old[j]);
      if (strcmp(penalty, "geLasso")==0) sG = sG + fabs(beta[l*p+j]) - fabs(beta_old[j]);
      if (strcmp(penalty, "geMCP")==0) sG = sG + MCP(beta[l*p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
      if (strcmp(penalty, "gMCP")==0) sG = sG + MCP(beta[l*p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
    }

    // Update df
    df[0] = df[0] + fabs(beta[l*p+j]) / fabs(z);
    }
}

// Groupwise local coordinate descent updates -- binomial
static void gLCD_binomial(double *beta, char *penalty, double *x, double *r, int g, int *K1, int n, int l, int p, double lam1, double lam2, double gamma, double tau, double *df, double *beta_old, double *w, double delta)
{
  // Calculate v
  int K = K1[g+1] - K1[g];
  double *v = Calloc(K, double);
  for (int j=K1[g]; j<K1[g+1]; j++) {
    for (int i=0; i<n; i++) v[j-K1[g]] += x[j*n+i] * w[i] * x[j*n+i];
    v[j-K1[g]] = v[j-K1[g]]/n;
  }

  // Make initial local approximation
  double sG = 0; // Sum of inner penalties for group
  if (strcmp(penalty, "geLasso")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(beta_old[j]);
  if (strcmp(penalty, "geMCP")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(beta_old[j], lam1, gamma);
  if (strcmp(penalty, "gMCP")==0) for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + MCP(beta_old[j], lam1, gamma);
  if (strcmp(penalty, "gBridge")==0) {
      for (int j=K1[g]; j<K1[g+1]; j++) sG = sG + fabs(beta_old[j]);
      if (sG==0) return;
      if (sG < delta) {
	for (int j=K1[g]; j<K1[g+1]; j++) {
	  beta[l*p+j] = 0;
	  for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
	}
	return;
    }
  }

  // Coordinate descent
  for (int j=K1[g]; j<K1[g+1]; j++) {

    // Calculate u
    double u = 0;
    for (int i=0; i<n; i++) u = u + x[j*n+i] * w[i] * r[i];
    u = u/n + v[j-K1[g]]*beta_old[j];

    // Calculate ljk
    double ljk=0;
    if (lam1 != 0) {
      if (strcmp(penalty, "gMCP")==0) ljk = dMCP(sG, lam1, (K*gamma*pow(lam1,2))/(2*lam1)) * dMCP(beta[l*p+j], lam1, gamma);
      if (strcmp(penalty, "geLasso")==0) ljk = lam1*exp(-tau/lam1*sG);
      if (strcmp(penalty, "geMCP")==0) ljk = exp(-tau/pow(lam1,2)*sG)*dMCP(beta[j],lam1,gamma);
      if (strcmp(penalty,"gBridge")==0) ljk = lam1 * gamma * pow(sG, gamma-1);
    }

    // Update beta
    beta[l*p+j] = S(u, ljk) / (v[j-K1[g]]*(1+lam2));
    //Rprintf("u: %f  ljk:%f  b:%f  oldbeta: %f\n", u, ljk, beta[l*p+j], beta_old[j]);

    // Update r
    if (beta[l*p+j] != beta_old[j]) {
      for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[n*j+i];
      if (strcmp(penalty, "gBridge")==0) sG = sG + fabs(beta[l*p+j]) - fabs(beta_old[j]);
      if (strcmp(penalty, "geLasso")==0) sG = sG + fabs(beta[l*p+j]) - fabs(beta_old[j]);
      if (strcmp(penalty, "geMCP")==0) sG = sG + MCP(beta[l*p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
      if (strcmp(penalty, "gMCP")==0) sG = sG + MCP(beta[l*p+j], lam1, gamma) - MCP(beta_old[j], lam1, gamma);
    }

    // Update df
    df[0] = df[0] + fabs(beta[l*p+j]) / fabs(u/v[j-K1[g]]);
    }
  Free(v);
}

static void gpPathFit_gaussian(double *beta, int *iter, double *df, double *loss, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, int *dfmax_, double *group_multiplier, int *user_)
{
  // Initialization of variables
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double delta=delta_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; int user=user_[0];
  double *r = Calloc(n, double);
  double *beta_old = Calloc(p, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=0; j<p; j++) {
      beta_old[j] = crossprod(x, r, n, j)/n;
      for (int i=0; i<n; i++) r[i] = r[i] - beta_old[j] * x[j*n+i];
    }
  }

  // Path
  for (int l=0; l<L; l++) {
    if (l != 0) for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l]++;

      // Check dfmax
      int a = 0;
      for (int j=0; j<p; j++) if (beta[l*p+j]!=0) a++;
      if (a > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	return;
      }
      df[l] = 0;

      // Update unpenalized covariates
      for (int j=0; j<K0; j++) {
	double shift = crossprod(x, r, n, j)/n;
	beta[l*p+j] = shift + beta_old[j];
	for (int i=0; i<n; i++) r[i] -= shift * x[n*j+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	if (user | l!=0) gLCD_gaussian(beta, penalty, x, r, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, &df[l], beta_old, delta);
      }

      // Check for convergence      
      if (checkConvergence(beta, beta_old, eps, l, p)) {
	converged  = 1;
	loss[l] = gLoss(r,n);
	break;
       }
      for (int j=0; j<p; j++) beta_old[j] = beta[l*p+j];
    }
  }
  Free(beta_old);
  Free(r);
}

static void gpPathFit_binomial(double *beta0, double *beta, int *iter, double *df, double *Dev, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, double *tau_, double *group_multiplier, int *dfmax_, int *warn_, int *user_)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double delta=delta_[0]; double gamma=gamma_[0]; double tau=tau_[0]; int dfmax=dfmax_[0]; int warn = warn_[0]; int user = user_[0];
  double beta0_old;
  double *r = Calloc(n, double);
  double *w = Calloc(n, double);
  double *beta_old = Calloc(p, double);

  // Initialization
  double ybar = mean(y,n);
  double nullDev = 0;
  for (int i=0;i<n;i++) nullDev += - y[i]*log(ybar) - (1-y[i])*log(1-ybar);
  if (strcmp(penalty, "gBridge")==0) {
    for (int j=0; j<p; j++) {
      double z=0;
      for (int i=0; i<n; i++) z = z + 0.25 * x[j*n+i] * (y[i]-0.5);
      beta_old[j] = z/n;
      for (int i=0; i<n; i++) r[i] = r[i] - beta_old[j] * x[j*n+i];
    }
  }

  // Path
  double eta, pi;
  for (int l=0; l<L; l++) {
    if (l != 0) {
      beta0_old = beta0[l-1];
      for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    }
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l] = iter[l] + 1;

      // Check dfmax
      int a = 0;
      for (int j=0; j<p; j++) if (beta[l*p+j]!=0) a++;
      if (a > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	return;
      }

      // Approximate L
      Dev[l] = 0;
      for (int i=0; i<n; i++) {
	eta = beta0_old;
	for (int j=0; j<p; j++) if (beta_old[j] != 0) eta = eta + x[j*n+i] * beta_old[j];
	pi = exp(eta) / (1+exp(eta));
	if (strncmp(penalty, "gr", 2)==0) r[i] = (y[i] - pi) / 0.25;
	else {
	  if (pi > .9999) {
	    pi = 1;
	    w[i] = .0001;
	  } else if (pi < .0001) {
	    pi = 0;
	    w[i] = .0001;
	  } else w[i] = pi*(1-pi);
	  r[i] = (y[i] - pi) / w[i];
	}
	Dev[l] = Dev[l] - y[i]*log(pi) - (1-y[i])*log(1-pi);
      }

      // Check for saturation
      if (Dev[l]/nullDev < .01) {
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

      // Update intercept
      double xwr = 0;
      double xwx = 0;
      for (int i=0;i<n;i++) {
	xwr = xwr + w[i]*r[i];
	xwx = xwx + w[i];
      }
      beta0[l] = xwr/xwx + beta0_old;
      for (int i=0;i<n;i++) r[i] = r[i] - (beta0[l] - beta0_old);
      df[l] = 1;
  
      // Update unpenalized covariates
      for (int j=0; j<K0; j++) {
	xwr = 0;
	xwx = 0;
	for (int i=0; i<n; i++) {
	  xwr = xwr + x[j*n+i]*w[i]*r[i];
	  xwx = xwx + x[j*n+i]*w[i]*x[j*n+i];
	}
	beta[l*p+j] = xwr/xwx + beta_old[j];
        for (int i=0; i<n; i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j]) * x[j*n+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	if (user | l!=0) gLCD_binomial(beta, penalty, x, r, g, K1, n, l, p, lam1[l]*group_multiplier[g], lam2[l], gamma, tau, &df[l], beta_old, w, delta);
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

static void grPathFit_gaussian(double *beta, int *iter, double *df, double *loss, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *group_multiplier, int *dfmax_, int *user_)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double gamma=gamma_[0]; int dfmax=dfmax_[0]; int user = user_[0];
  double *r = Calloc(n, double);
  double *beta_old = Calloc(p, double);
  int *active = Calloc(J, int);

  // Initialization
  for (int i=0; i<n; i++) r[i] = y[i];

  // Path
  for (int l=0; l<L; l++) {
    if (l != 0) for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l]++;

      // Check dfmax
      int a = 0;
      for (int j=0; j<J; j++) a += (beta_old[j] !=0);
      if (a > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	Free(active);
	return;
      }
      df[l] = 0;

      // Update unpenalized covariates
      for (int j=0; j<K0; j++) {
	double shift = crossprod(x, r, n, j)/n;
	beta[l*p+j] = shift + beta_old[j];
	for (int i=0; i<n; i++) r[i] -= shift * x[n*j+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	//Rprintf("g=%d  K1[g]=%d  K1[g+1]=%d  gm[g] = %f\n",g,K1[g],K1[g+1],group_multiplier[g]);
	if (user | l!=0) gd_gaussian(beta, x, r, g, K1, active, n, l, p, penalty, lam1[l]*group_multiplier[g], lam2[l], gamma, &df[l], beta_old);
      }

      // Check convergence
      if (checkConvergence(beta, beta_old, eps, l, p)) {
	converged  = 1;
	loss[l] = gLoss(r,n);
	break;
      }
      for (int j=0; j<p; j++) beta_old[j] = beta[l*p+j];
    }
  }
  Free(beta_old);
  Free(active);
  Free(r);
}

static void grPathFit_binomial(double *beta0, double *beta, int *iter, double *df, double *Dev, double *x, double *y, int *n_, int *p_, char **penalty_, int *J_, int *K1, int *K0_, double *lam1, double *lam2, int *L_, double *eps_, int *max_iter_, double *gamma_, double *group_multiplier, int *dfmax_, int *warn_, int *user_)
{
  int n=n_[0]; int p=p_[0]; char *penalty=penalty_[0]; int J=J_[0]; int K0=K0_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double gamma=gamma_[0]; int dfmax=dfmax_[0]; int warn = warn_[0]; int user = user_[0];
  double beta0_old;
  double *r = Calloc(n, double);
  double *beta_old = Calloc(p, double);
  int *active = Calloc(J, int);

  // Initialization
  double ybar = mean(y,n);
  double nullDev = 0;
  for (int i=0; i<n; i++) nullDev = nullDev - y[i]*log(ybar) - (1-y[i])*log(1-ybar);

  // Path
  double eta, pi;
  for (int l=0; l<L; l++) {
    if (l != 0) {
      beta0_old = beta0[l-1];
      for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
    }
    while (iter[l] < max_iter) {
      int converged = 0;
      iter[l]++;

      // Check dfmax
      int a = 0;
      for (int j=0; j<J; j++) a += (beta_old[j] !=0);
      if (a > dfmax) {
	for (int ll=l; ll<L; ll++) {
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(r);
	Free(active);
	return;
      }

      // Approximate L
      Dev[l] = 0;
      for (int i=0; i<n; i++) {
	eta = beta0_old;
	for (int j=0; j<K0; j++) eta += x[n*j+i] * beta_old[j];
	for (int g=0; g<J; g++) {
	  if (active[g]) {
	    for (int j=K1[g]; j<K1[g+1]; j++) eta += x[n*j+i] * beta_old[j];
	  }
	}
	pi = exp(eta) / (1+exp(eta));
	r[i] = (y[i] - pi) / 0.25;
	Dev[l] += - y[i]*log(pi) - (1-y[i])*log(1-pi);
	}

      // Check for saturation
      if (Dev[l]/nullDev < .01) {
	if (warn) warning("Model saturated; exiting...");
	for (int ll=l; ll<L; ll++) {
	  beta0[ll] = R_NaReal;
	  for (int j=0; j<p; j++) beta[ll*p+j] = R_NaReal;
	}
	Free(beta_old);
	Free(active);
	Free(r);
	return;
      }

      // Update intercept
      double shift = mean(r,n);
      beta0[l] = shift + beta0_old;
      for (int i=0; i<n; i++) r[i] -= shift;
      df[l] = 1;

      // Update unpenalized covariates
      for (int j=0; j<K0; j++) {
	shift = crossprod(x, r, n, j)/n;
	beta[l*p+j] = shift + beta_old[j];
	for (int i=0; i<n; i++) r[i] -= shift * x[n*j+i];
	df[l] = df[l] + 1;
      }

      // Update penalized groups
      for (int g=0; g<J; g++) {
	if (user | l!=0) gd_binomial(beta, x, r, g, K1, active, n, l, p, penalty, lam1[l]*group_multiplier[g], lam2[l], gamma, &df[l], beta_old);
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
  Free(active);
  Free(r);
}

static const R_CMethodDef cMethods[] = {
  {"gpPathFit_gaussian", (DL_FUNC) &gpPathFit_gaussian, 23},
  {"gpPathFit_binomial", (DL_FUNC) &gpPathFit_binomial, 25},
  {"grPathFit_gaussian", (DL_FUNC) &grPathFit_gaussian, 21},
  {"grPathFit_binomial", (DL_FUNC) &grPathFit_binomial, 23},
  NULL
};

void R_init_grpreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
