#include <math.h>
#include <R.h>
#include <R_ext/Random.h>
#include "QRM.h"

void frank(int *n, double *theta, double *output)
{
  int i, k;
  double U, alpha, thetaval, p;
  thetaval = *theta;
  alpha = 1-exp(-thetaval);
  RANDIN; 

  for( i=0; i<*n; i++ ) {
    k = 1;
    U = unif_rand();
    p = alpha/thetaval;
    while(U > p) {
      k++;
      p = p + pow(alpha,k)/(k*thetaval);
    }
    *(output+i)=k;
  }
  RANDOUT;
}

void rgig(int *n, double *r, double *s, double *p, double *k1, double *k2, double *lambda, double *chi, double *psi, double *s1, double *s2, double *xsim)
{
  int i=0;
  int count =0;
  double U, Ustar, level, x;
  RANDIN;

  while (i < *n){
    U = unif_rand();
    Ustar = unif_rand();

    count++;
    if (U <= *r){
      x = log(1+(*s)*U/(*k1))/(*s);  
      level = log(ef(x, *lambda, *chi,*psi+2*(*s))/(*s1));
      if (log(Ustar) <= level){
	*(xsim+i) = x;
	i++;
      }
    } else {
      x = -log((*p)*(1-U)/(*k2))/(*p);
      level = log(ef(x, *lambda, *chi,*psi-2*(*p))/(*s2));
      if (log(Ustar) <= level){
	*(xsim+i) = x;
	i++;
      }
    }

  }
  *n = count;
  RANDOUT;
}

double ef(double x, double lambda, double chi, double psi)
{
  double result;
  result = pow(x,lambda-1.0)*exp(-0.5*(psi*x + chi/x));
  return result;
}

void SEprocExciteFunc(int*n, double *times, int *nmarks, double *marktimes, double *marks, double *beta, int *model, double *result)
{
  int i=0, j;
  double thetime, gamma, delta, rho=0.0, tmp;

  gamma = *beta;
  delta = 0.0;
  if (*model ==2 )
    /* Hawkes with mark influence */
    delta = *(beta+1);
  if (*model == 3)
    /* ETAS without mark influence */
    rho = *(beta+1);
  if (*model == 4){
    /* ETAS with mark influence */
    rho = *(beta+1);
    delta = *(beta+2);
  }
  while (i < *n){
    tmp = 0.0;
    thetime = *(times+i);
    j = 0;
    while ((*(marktimes+j) < thetime) & (j < *nmarks)){
      if (*model == 1)
	tmp = tmp + contribH((thetime-*(marktimes+j)),0.0,gamma,delta);
      if (*model == 2)
	tmp = tmp + contribH((thetime-*(marktimes+j)),*(marks+j),gamma,delta);
      if (*model == 3)
	tmp = tmp + contribE((thetime-*(marktimes+j)),0.0,gamma,rho,delta);
      if (*model == 4)
	tmp = tmp + contribE((thetime-*(marktimes+j)),*(marks+j),gamma,rho,delta);
      j++;
    }
    *(result+i) = tmp;
    i++;
  }
}

double contribH(double s, double y, double gammaval, double deltaval)
{
  double result;
  result = (1+deltaval*y)*exp(-gammaval*s);
  return result;
}

double contribE(double s, double y, double gammaval, double rhoval, double deltaval)
{
  double result;
  result = (1+deltaval*y)/pow(1+s/gammaval,(1.0+rhoval));
  return result;
}



