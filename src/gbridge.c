#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n);
static void free_vector(double *v);
static double **matrix(int nr, int nc);
static void free_matrix(double **M, int nr);
static double **as_matrix(double *v, int nr, int nc);
static void as_vector(double *v, double **M, int nr, int nc);

static double *vector(int n) {
  double *v;
  v = Calloc(n, double);
  return v;
}

static void free_vector(double *v) {
  Free(v);
}

static double **matrix(int nr, int nc) {
  int   i;
  double **M;
  M = Calloc(nr, double *);
  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  return M;
}

static void free_matrix(double **M, int nr) {
  int   i;
  for (i = 0; i < nr; i++) Free(M[i]);
  Free(M);
}

static double **as_matrix(double *v, int nr, int nc) {
  int i,j;
  double **M;

  M = Calloc(nr, double *);

  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) M[i][j] = v[j*nr+i];
  }
  return M;
}

static double **t_as_matrix(double *v, int nr, int nc)
{
  int i,j;
  double **M;

  M = Calloc(nr, double *);

  for (i = 0; i < nr; i++) M[i] = Calloc(nc, double);
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) M[i][j] = v[j*nr+i];
  }
  return M;
}

static void as_vector(double *v,double **M, int nr, int nc) {
  int i,j;

  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) v[j*nr+i] = M[i][j];
  }
}

static int checkConvergence(double *beta_new, double *beta_old, double eps, int p) {
  int j;
  int converged = 1;
  for (j=0; j < p; j++) {
    if (beta_new[j]!=0 & beta_old[j]!=0) {
      if (fabs((beta_new[j]-beta_old[j])/beta_old[j]) > eps) {
	converged = 0;
	break;
      } 
    } else if (beta_new[j]==0 & beta_old[j]!=0) {
      converged = 0;
      break;
    } else if (beta_new[j]!=0 & beta_old[j]==0) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

static double S(double x,double y)
{
  if (fabs(x) <= y) return(0); else {
    if (x > 0) return(x-y);
    else return(x+y);
  }
}

static void gLasso(double *beta, double *x, double *w, double *r, int K0, int Kj, int n, double lambda, double delta, double lambda2, double *df)
{
  int i, j, k, j1, j2, k1, k2, K;
  K = Kj - K0;
  double sxr, sxx, oldbeta, gradient_norm, sbb, ljk, s;
  double *u;
  u = vector(K);

  for (j=K0; j<Kj; j++) sbb = sbb + pow(beta[j],2);
  if (sbb==0)
    {
      gradient_norm = 0;
      for (j1=K0; j1<Kj; j1++)
	{
	  u[j1-K0] = 0;
	  for (i=0; i<n; i++) u[j1-K0] = u[j1-K0] + x[n*j1+i]*w[i]*r[i];
	  gradient_norm = gradient_norm + pow(u[j1-K0],2);
	}
      gradient_norm = sqrt(gradient_norm);
    }
  else
    {
      gradient_norm = 0;
      for (j1=K0; j1<Kj; j1++)
	{
	  u[j1-K0] = 0;
	  for (i=0; i<n; i++)
	    {
	      u[j1-K0] = u[j1-K0] + x[n*j1+i]*w[i]*r[i];
	      for (j2=K0; j2<Kj; j2++)
		{
		  u[j1-K0] = u[j1-K0] + x[n*j1+i]*w[i]*x[n*j2+i]*beta[j2];
		}
	    }
	  gradient_norm = gradient_norm + pow(u[j1-K0],2);
	}
      gradient_norm = sqrt(gradient_norm);
    }
  /*if (gradient_norm/n > sqrt(K)*lambda)*/
  if (gradient_norm/n > lambda)
    {
      sbb = sbb + delta;
      for (j=K0; j<Kj; j++)
	{
	  oldbeta = beta[j];
	  sxr=0;
	  for (i=0; i<n; i++) sxr = sxr + w[i]*x[n*j+i]*r[i];
	  if (w[0]==1) sxx=n;
	  else
	    {
	      sxx=0;
	      for (i=0; i<n; i++) sxx = sxx + w[i]*pow(x[n*j+i],2);
	    }
	  ljk = lambda/sqrt(sbb);
	  beta[j] = (sxr+sxx*oldbeta)/(sxx+n*(ljk + lambda2));
	  for (i=0; i<n; i++) r[i] = r[i] - (beta[j]-oldbeta)*x[n*j+i];
	  sbb = sbb + pow(beta[j],2) - pow(oldbeta,2);
	  df[0] = df[0] + fabs(beta[j]) / fabs(sxr/sxx+beta[j]);
	}
    }
  else if (beta[K0]!=0)
    {
      for (j=K0; j<Kj; j++)
	{
	  oldbeta = beta[j];
	  beta[j] = 0;
	  for (i=0; i<n; i++) r[i] = r[i] + oldbeta*x[n*j+i];
	}
    }
  free_vector(u);
}

static double MCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta <= a*l) return(l*theta - pow(theta,2)/(2*a));
  else return(a*pow(l,2)/2);
}

static double dMCP(double theta, double l, double a)
{
  theta = fabs(theta);
  if (theta < a*l) return(l-theta/a);
  else return(0);
}

static void gLCD(double *beta, char *penalty, double *x, double *w, double *r, int K0, int Kj, int n, double lambda, double lambda2, double gamma, double tau, double *df)
{
  int K = Kj - K0;
  double sG = 0;
  double delta = 0.0000001;
  double z, v, oldbeta, ljk;

  if (strcmp(penalty,"gBridge")==0 | strcmp(penalty,"geLasso")==0)
    {
      for (int j=K0; j<Kj; j++) sG = sG + fabs(beta[j]);
      if (sG==0 & strcmp(penalty,"gBridge")==0) return;
    }
  if (strcmp(penalty,"gMCP")==0 | strcmp(penalty,"geMCP")==0) for (int j=K0; j<Kj; j++) sG = sG + MCP(beta[j],lambda,gamma);

  for (int j=K0; j<Kj; j++)
    {
      if (sG < delta & strcmp(penalty,"gBridge")==0) beta[j] = 0;
      else
	{
	  oldbeta = beta[j];
	  z=v=0;
	  for (int i=0; i<n; i++) z = z + w[i]*x[n*j+i]*r[i];
	  if (w[0]==1) v=1;
	  else
	    {
	      for (int i=0; i<n; i++) v = v + w[i]*pow(x[n*j+i],2);
	      v = v/n;
	    }
	  z = z/n + v*oldbeta;
	  if (lambda==0) ljk=0;
	  else
	    {
	      if (strcmp(penalty,"gBridge")==0) ljk = sqrt(K)*lambda*gamma*pow(sG,gamma-1);
	      if (strcmp(penalty,"gMCP")==0) ljk = dMCP(sG,lambda,(K*gamma*pow(lambda,2))/(2*lambda))*dMCP(beta[j],lambda,gamma);
	      if (strcmp(penalty,"geLasso")==0) ljk = lambda*exp(-tau/lambda*sG);
	      if (strcmp(penalty,"geMCP")==0) ljk = exp(-tau/pow(lambda,2)*sG)*dMCP(beta[j],lambda,gamma);
	    }
	  beta[j] = S(z,ljk)/(v+lambda2);
	}
      if (beta[j]!=oldbeta)
	{
	  for (int i=0; i<n; i++) r[i] = r[i] - (beta[j]-oldbeta)*x[n*j+i];
	  if (strcmp(penalty,"gBridge")==0 | strcmp(penalty,"geLasso")==0) sG = sG + fabs(beta[j]) - fabs(oldbeta);
	  if (strcmp(penalty,"gMCP")==0 | strcmp(penalty,"geMCP")==0) sG = sG + MCP(beta[j],lambda,gamma) - MCP(oldbeta,lambda,gamma);
	}
      df[0] = df[0] + fabs(beta[j])/fabs(z/v);
    }
}

static void gpPathFit(double *beta, int *counter, double *df, double *x, double *y, int *group, char **family, int *n, int *p, int *J, int *K, char **penalty, double *lambda, int *nlambda, double *eps, int *max_iter, int *verbose, double *delta, double *gamma, double *tau, double *lambda2, double *group_multiplier, int *dfmax, int *warn_conv)
{
  int l, i, j, g, K0, Kj, converged, active;
  int saturated=0;
  double sxr, sxx, oldbeta;
  double ybar, yp, yy;
  double eta, pi;
  double **Beta, *r, *w, *beta_old;
  Beta = matrix(*nlambda,*p);
  r = vector(*n);
  beta_old = vector(*p);
  w = vector(*n);

  /* Initial setup */
  if (strcmp(penalty[0],"gBridge")==0) {
    if (strcmp(family[0],"gaussian")==0) {
      if (group[0]==0) {
	sxr=0;
	for (i=0; i<*n; i++) sxr = sxr + y[i];
	beta_old[0] = sxr / *n;
	for (i=0;i<*n;i++) {
	  w[i] = 1;
	  r[i] = y[i] - beta_old[0];
	}
	j=1;
      } else for (i=0;i<*n;i++) {
	  w[i] = 1;
	  r[i] = y[i];
	  j=0;
	}
      for (;j<*p;j++) {
	sxr=0;
	for (i=0; i<*n; i++) sxr = sxr + x[j*n[0]+i]*r[i];
	beta_old[j] = sxr / *n;
      }
      for (i=0;i<n[0];i++) {
	r[i] = y[i];
	for (j=0;j<*p;j++) r[i] = r[i] - x[j*n[0] + i]*beta_old[j];
      }
    } else {
      if (group[0]==0) {
	sxr=0;
	for (i=0; i<*n; i++) sxr = sxr + y[i];
	pi = sxr/n[0];
	beta_old[0] = log((pi)/(1-pi));
	for (i=0;i<*n;i++) {
	  w[i] = sqrt(pi*(1-pi));
	  r[i] = (y[i] - pi)/w[i];
	}
	j=1;
      } else for (i=0;i<*n;i++) {
	  w[i] = 1;
	  r[i] = (y[i] - 0.5)/0.5;
	  j=0;
	}
      for (;j<*p;j++) {
	sxr=0;
	for (i=0; i<*n; i++) sxr = sxr + w[i]*x[j*n[0]+i]*r[i];
	beta_old[j] = sxr / *n;
      }
    }
  } else {
    for (j=0;j<*p;j++) beta_old[j] = 0;
    if (strcmp(family[0],"gaussian")==0) {
      for (i=0;i<*n;i++) {
	w[i] = 1;
	r[i] = y[i];
      }
    }
  }

  /* Path */
  for (l=0;l<*nlambda;l++)
    {
      if (l==0) for (j=0;j<*p;j++) Beta[0][j] = beta_old[j];
      else for (j=0;j<*p;j++) Beta[l][j] = beta_old[j] = Beta[l-1][j];

      /* Check dfmax */
      active = 0;
      for (int j=0;j<*p;j++) if (Beta[l][j]!=0) active++;
      if (active > *dfmax)
	{
	  for (int ll=l;ll<*nlambda;ll++)
	    {
	      for (int j=0;j<*p;j++) Beta[ll][j] = R_NaReal;
	    }
	  as_vector(beta,Beta,*nlambda,*p);
	  free_matrix(Beta,*nlambda);
	  free_vector(beta_old);
	  free_vector(w);
	  free_vector(r);
	  return;
	}

      if (*verbose) Rprintf("Starting new fit: lambda = %f\n",lambda[l]);
      while (counter[l] < *max_iter)
	{
	  converged = 0;
	  counter[l] = counter[l] + 1;
	  if (*verbose) Rprintf("Iteration: %d\n",counter[l]);

	  /* Approximate L                 */
	  if (strcmp(family[0],"binomial")==0)
	    {
	      ybar = 0;
	      yp = 0;
	      yy = 0; for (i=0;i<*n;i++) ybar = ybar + y[i];
	      ybar = ybar / *n;
	      for (i=0;i<*n;i++)
		{
		  eta = 0;
		  for (j=0;j<*p;j++) eta = eta + x[j * *n + i]*Beta[l][j];
		  pi = exp(eta)/(1+exp(eta));
		  if (pi > .999)
		    {
		      pi = 1;
		      w[i] = .001;
		    }
		  else if (pi < .001)
		    {
		      pi = 0;
		      w[i] = .001;
		    }
		  else w[i] = sqrt(pi*(1-pi));
		  r[i] = (y[i] - pi)/w[i];
		  yp = yp + pow(y[i]-pi,2);
		  yy = yy + pow(y[i]-ybar,2);
		}
	    }
	  if (strcmp(family[0],"binomial")==0)
	    {
	      if (yp/yy < .01)
		{
		  warning("Model saturated; exiting...");
		  saturated=1;
		  break;
		}
	    }
	  df[l] = 0;
  
	  /* Update unpenalized covariates */
	  for (j=0;; j++)
	    {
	      if (group[j]!=0) break;
	      sxr=0;
	      for (i=0; i<*n; i++) sxr = sxr + w[i]*x[j * *n + i]*r[i];
	      if (w[0]==1) sxx=*n;
	      else
		{
		  sxx=0;
		  for (i=0; i<*n; i++) sxx = sxx + w[i]*pow(x[j * *n + i],2);
		}

	      oldbeta = Beta[l][j];
	      Beta[l][j] = sxr/sxx + Beta[l][j];
	      for (i=0; i<*n; i++) r[i] = r[i] - (Beta[l][j]-oldbeta)*x[j * *n + i];
	      df[l] = df[l] + 1;
	    }
	  /* Update penalized covariates */
	  for (g=0; g<*J; g++)
	    {
	      K0 = j;
	      Kj = j + K[g];
	      if (strcmp(penalty[0],"gLasso")==0) gLasso(Beta[l],x,w,r,K0,Kj,*n,lambda[l]*group_multiplier[g],*delta,lambda2[l],&df[l]);
	      else gLCD(Beta[l],penalty[0],x,w,r,K0,Kj,*n,lambda[l]*group_multiplier[g],lambda2[l],*gamma,*tau,&df[l]);
	      j = Kj;
	    }
	  if (checkConvergence(Beta[l],beta_old,*eps,*p))
	    {
	      converged  = 1;
	      break;
	    }
	  for (j=0;j<*p;j++) beta_old[j] = Beta[l][j];
	}
      if (saturated) break;
      if (converged==0 & warn_conv[0]) warning("Failed to converge");
    }
  as_vector(beta,Beta,*nlambda,*p);

  free_matrix(Beta,*nlambda);
  free_vector(beta_old);
  free_vector(w);
  free_vector(r);
}

static const R_CMethodDef cMethods[] = {
  {"gpPathFit", (DL_FUNC) &gpPathFit, 24},
  NULL
};

void R_init_grpreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
