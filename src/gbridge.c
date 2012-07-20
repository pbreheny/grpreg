#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n);
static int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
static double gLoss(double *r, int n);
static double S(double z, double l);

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

static void gBridgeFit_gaussian(double *beta, int *iter, double *df, double *loss, double *x, double *y, int *group, int *n_, int *p_, int *J_, int *K, double *lam1, double *lam2, int *L_, double *eps_, double *delta_, int *max_iter_, double *gamma_, int *dfmax_, int *group_multiplier)
{
  // Initialization of variables
  int n=n_[0]; int p=p_[0]; int J=J_[0]; int L=L_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double delta=delta[0]; double gamma=gamma_[0]; int dfmax=dfmax_[0]; 
  double *r, *beta_old;
  r = vector(n);
  for (int i=0; i<n; i++) r[i] = y[i];
  beta_old = vector(p);

  // Path
  for (int l=0; l<L; l++) {
    //Rprintf("lambda: %f %f\n",lam1[l], lam2[l]);
    if (l != 0) for (int j=0; j<p; j++) beta_old[j] = beta[(l-1)*p+j];
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
      df[l] = 0;

      // Update unpenalized covariates
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
	loss[l] = gLoss(r,n);
	break;
       }
      for (j=0; j<p; j++) beta_old[j] = beta[l*p+j];
    }
  }
  Free(beta_old);
  Free(r);
}

static const R_CMethodDef cMethods[] = {
  {"gpPathFit", (DL_FUNC) &gpPathFit, 24},
  NULL
};

void R_init_grpreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
