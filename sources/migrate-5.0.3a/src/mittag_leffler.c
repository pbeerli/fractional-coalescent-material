/*-----------------------------------------------------------------
  Bayesian inference of population genetic forces: drift, migraiton, divergence
  allowing for the n-coalescent, the f-coalescent, and the BSC-coalescent
 
  Peter Beerli
  Department of Scientific Computing
  Florida State University
  Tallahassee FL 32306-4120
  beerli@fsu.edu
 
  Copyright 2017 Peter Beerli, Tallahassee FL

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*-----------------------------------------------------------------
*/

// to test standalone use
//gcc -g -DSTANDALONEMITTAGLEFFLER mittag_leffler.c hermite_interpoly.c romberg.c mittag_leffler_interpol_data.c -o mlf

#include "mittag_leffler.h"
#include <math.h>
#include "migration.h"
#include "random.h"
#include "romberg.h"
#include "hermite_interpoly.h"
#include "sighandler.h"
#include "speciate.h"

//!Mittage-Leffler function in the real case
//!Somayeh Mashayekhi March 2017
// translated from Fortran to C Peter Beerli March 2017
#include <complex.h>    // Standard Library of Complex Numbers
//#undef I
#define J ((complex double) I)  
#ifndef Pi
#define Pi 3.14159265358979323846264338328
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

complex double mittag_leffler(double alpha, double beta, complex double z);
complex double  K(double alpha, double beta, double x, complex double z) ;
complex double P(double alpha, double beta, double eps, double phi, complex double z); // returns , complex double P1
complex double  ML(complex double z, double alpha, double beta, double Q, double X0);
double KK(double x, complex double *args, long arglen);
double PP(double phi, complex double *args, long arglen);
double interval_mittag_leffler(double r, double alpha, double lambda, double tmin, double tmax);
void set_mittag_leffler(option_fmt * options);
double interval_mittag_leffler_func(double r, double alpha, double t0, double mu, double sigma, species_fmt *s, double tmin, double tmax);
double propose_new_mlftime(double lambda, double alpha, double r1, double r2);


// calculates internals of the generalized mittag-leffler function for z with alpha and beta
// with precision Q and X0
complex double  ML(complex double z, double alpha, double beta, double Q, double X0) 
{
  long K0;
  complex double b = 0.0;
  complex double kargs[3] = {alpha, beta, z};
  complex double pargs[4] = {alpha, beta, 1.0, z};

  if (1.0 < alpha) //alpha > 1;  beta in Real;  z in Complex; 
    {
      K0 = (long) (floor(alpha) + 1.0);
      int i;
      for (i=0; i<K0; ++i)
	{
	  b += cpow(z,1.0/K0) * cexp(2.0*Pi*J*i/K0);
	}
      //printf("1 z= %.5f alpha > 1.0 \n", creal(z));
      return b;
    }
  else if (alpha == 1.0  && beta == 1.0)
    {
      b = cexp(z);
      return b;
    }
  else if (z==0.0) // z==0
    {
      //printf("2 (z=%.5f)==0 alpha < 1.0 \n", creal(z));
      b = 1.0/tgamma(beta);
      return b;
    }
  else if (cabs(z)<1.0) // |z| < 1.0 z:{-1,1}
    {
      K0 = (long) max(creal((ceil ((1.0-beta)/alpha))),ceil(creal(clog(Q*(1.0-cabs(z)))/clog(cabs(z)))));
      K0 = min(200,K0);
      b = 0.0;
      int i;
      for (i=0;i<K0+1;i++)
	{
	  b += (cpow(z,i))/ tgamma(beta+alpha*i);
	}
      return b;
    }
  else if (cabs(z)>floor(10.0+5.0*alpha))
    {
      K0= (long) (floor(creal(-clog(Q)/clog(cabs(z)))));
      K0 = min(200,K0);
      double c = 0.0;
      int i;
      for(i=1; i<K0+1; i++)
	{
	  c += cpow(z,(-i)) / tgamma(beta-alpha*i);
	}
      if(fabs(carg(z)) < (alpha*Pi/4.0 + 0.5 * min(Pi, alpha*Pi)))
	{
	  b=(1.0/alpha)*cpow(z,((1-beta)/alpha))*cexp(cpow(z,(1.0/alpha))) - c;
	}
      else 
	{
	  b = -c;
	}
      return b;
    }
  else
    {
      // compute X0 (done outside)
      if(fabs(carg(z)) > alpha * Pi)
	{
	  if (beta <= 1.0)
	    {
	      b = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
	    }
	  else
	    {
	      b = romberg(&KK, kargs, 3, 1.0, X0, 100, Q);
	      b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
	    }
	  return b;
	}
      else if(fabs(carg(z)) < alpha * Pi)
	{
	  if (beta <= 1)
	    {
	      b = romberg(&KK, kargs, 3, 0.0, X0, 100, Q);
	      b += (1/alpha)*cpow(z,((1.0-beta)/alpha)) * cexp(cpow(z,(1.0/alpha)));
	    }
	  else
	    {
	      pargs[2]= cabs(z)/2.0;
	      b = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	      b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
	      b += (1/alpha)*cpow(z,((1-beta)/alpha))*cexp(cpow(z,(1/alpha)));
	    }
	  return b;
	}
      else
	{
	  double lowbound = creal((cabs(z) + 1.0)/ 2.0);
	  pargs[2]= lowbound;
	  b = romberg(&KK, kargs, 3, creal(pargs[2]) , X0, 100, Q);
	  b += romberg(&PP, pargs, 3, -alpha*Pi, alpha*Pi, 100, Q);
	}
      return b;
    }
}

complex double  K(double alpha, double beta, double x, complex double z) 
{
  if(x==0.0)
    return 0.0;
  complex double K1=(1/(alpha*Pi))*cpow(x,((1.0-beta)/alpha))*
    cexp(-cpow(x,(1.0/alpha)))*
    ((x*csin(Pi*(1.0-beta))-z*csin(Pi*(1.0-beta+alpha)))/(x*x-2.0*x*z*ccos(alpha*Pi)+z*z));
  return K1;
}

complex double P(double alpha, double beta, double eps, double phi, complex double z) // returns , complex double P1
{
  complex double omega = phi * (1.0+ (1.0-beta)/alpha) + cpow(eps,(1.0/alpha)) * csin(phi/alpha);    
  complex double P1 = (1.0/(2.0*alpha*Pi)) * cpow(eps,(1.0+(1.0-beta)/alpha)) * cexp(cpow(eps,(1.0/alpha)) * ccos(phi/alpha)) * (ccos(omega)+J*csin(omega))/(eps*cexp(J*phi)-z);
  return P1;
}

//double mlf(double alpha, double beta, complex double z)
complex double mittag_leffler(double alpha, double beta, complex double z)
{
  
  double Q=PRECISION;  // set in definitions.h
  double X0;

  if(beta >= 0.0)
    {
      X0 = max(
	       max(creal(1.0), creal(2.0*creal(z))),
	       creal(cpow((-clog(Pi*Q/creal(6.0))),(alpha))));
    }
  else
    {
      X0 = max(
	       max(creal(cpow((fabs(beta)+1.0),(alpha))), creal(2.0*creal(z))),
	       creal(cpow((-2.0*clog(Pi*Q/creal(6*(fabs(beta)+2.0)* cpow((2.0*fabs(beta)),fabs(beta))))),(alpha))));
    }
#ifndef MLF_SLOW
  complex double b;
  if (alpha == beta || beta==1.0)
    b = MLinterpol(z,alpha, beta, Q, X0);
  else
    b = clog(ML(z,alpha, beta, Q, X0));
#ifdef MLFCHECK
  complex double c = clog(ML(z,alpha, beta, Q, X0));
  printf("@CHECK@ interpol= %f mlf= %f z= %f a= %f b= %f\n",creal(b),creal(c),creal(z),alpha,beta);
#endif
#else
  complex double b = clog(ML(z,alpha, beta, Q, X0));
#endif
  return b;
}

double KK(double x, complex double *args, long arglen)
{
  (void) arglen;
  double alpha = (double) args[0];
  double beta  = (double) args[1];
  complex double z = args[2];
  complex double b = K(alpha,beta,x,z);
  return creal(b);
}

double PP(double phi, complex double *args, long arglen)
{
  (void) arglen;
  double alpha = (double) args[0];
  double beta  = (double) args[1];
  double eps   = (double) args[2];
  complex double z = args[3];
  complex double b = P(alpha,beta,eps,phi,z);
  return creal(b);
}
//1) For drawing time you need to solve
//   r=E_alpha(-lambda *t^alpha).
//   Before you had r=E_alpha(-lambda *t).

//2) For probability you need to write
//   t^(alpha-1)*E_alpha,alpha(-lambda *t^alpha).
//   Before you       had E_alpha,alpha(-lambda *t).


double interval_mittag_leffler(double r, double alpha, double lambda, double tmin, double tmax)
{
  double minu=tmin;
  double maxu=tmax;
  double u=(maxu+minu)/2.0;
  double prob;
  long maxcount=200;
  double beta = 1.0;
  double logr = log(r);
  while(maxu-minu>SMALLEPSILON && maxcount-- > 0)
    {
      prob =  mittag_leffler(alpha,beta, (complex double) (-(pow(u,alpha)) * lambda));
      //printf("@ %f %f %f prob=%f r=%f \n",minu,u,maxu, prob, r);
      if (logr>prob)
	{
	  maxu = u;
	}
      else
	{
	  minu = u;
	}
      u = (maxu+minu)/2.0;
    }
  if (maxcount<=0)
    {
      fprintf(stderr,"WARNING: drawing a random time for MLF failed: minu=%f u=%f maxu=%f\n",minu,u,maxu);
    }
  return u;
}

#ifndef STANDALONEMITTAGLEFFLER
double interval_mittag_leffler_func(double r, double alpha, double t0, double mu, double sigma, species_fmt *s, double tmin, double tmax)
{
  double minu=tmin;
  double maxu=tmax;
  double u=(maxu+minu)/2.0;
  double prob;
  long maxcount=200;
  double beta = 1.0;
  double logr = log(r);
  double lambda;
  while(maxu-minu>EPSILON && maxcount-- > 0)
    {
      double t1 = t0 + u;
      lambda = (*log_prob_wait_speciate)(t0,t1,mu,sigma,s);
      prob = mittag_leffler(alpha,beta, (complex double) (- pow(u,alpha) * lambda));
      //printf("@ %f %f %f prob=%f r=%f \n",minu,u,maxu, prob, r);
      if (logr>prob)
	{
	  maxu = u;
	}
      else
	{
	  minu = u;
	}
      u = (maxu+minu)/2.0;
    }
  if (maxcount<=0)
    warning("drawing a random time for MLF failed: minu=%f u=%f maxu=%f\n",minu,u,maxu);
  return u;
}

// Somayeh Mashayekhi 
double propose_new_mlftime(double lambda, double alpha, double r1, double r2)
{
  double pia = PI * alpha;
  double denoma = 1.0 / alpha;
  double denomlambda = 1.0 / lambda;
  //return -pow(denomlambda,denoma) * pow(sin(pia)/(tan(pia*(1-r1))) - cos(pia),denoma) * log(r2);
  return -pow(denomlambda * (sin(pia)/(tan(pia*(1-r1))) - cos(pia)),denoma) * log(r2);
}

void set_mittag_leffler(option_fmt * options)
{
  char *input;
  input = (char *) mycalloc(LINESIZE,sizeof(char));
  printf("Mittag-Leffler distribution\n");
  printf("----------------------------------------------------------------\n");
  printf("The standard Coalescent is based on the Exponential distribution\n");
  printf("Mittag-Leffler is a generalization with an additional parameter\n");
  printf("alpha, ranging between 0 and 1; With an alpha=1.0 we recover the\n");
  printf("standard Exponential; with alpha < 1 the distribution has\n");
  printf("shorter coalescences near today, values are allowed\n");
  printf("between 0.01 and 1.00, using 2-digit precision.\n");
  printf("[Default: alpha=1.00]\n");
  printf("> ");fflush(stdout);
  FGETS(input,LINESIZE,stdin);
  if(input[0]!='\0')
    {
      options->mlalpha = atof(input);
    }
  else
    {
      options->mlalpha = 1.0;
    }
}
#endif



#ifdef STANDALONEMITTAGLEFFLER
void run_grid(double alpha,double beta)
{
  complex double zmax = 100;
  complex double zmin = -100;
  complex double z;
  for (z=zmin; creal(z) < creal(zmax); z += 1.0)
    {
      complex double b = mittag_leffler(alpha, beta, z);
      printf("mlf(z=%f,alpha=%f,beta=%f) =  %.8g %+.8gi\n", creal(z), alpha, beta, creal(b), cimag(b));
    }
}


void help(int argc, char **argv);

void help(int argc, char **argv)
{
  if ((argc==1) || (argc>1 && argv[1][0]=='-' && argv[1][1]=='h'))
    {
      printf("Syntax: mlf g           #prints mlf(lambda=[-100..99],alpha=0.6,beta=0.6)\n");
      printf("Syntax: mlf z           #prints mlf(lambda=z,alpha=0.6,beta=0.6)\n");
      printf("Syntax: mlf z a         #prints mlf(lambda=z,alpha=a,beta=a)\n");
      printf("Syntax: mlf z a b       #prints mlf(lambda=z,alpha=a,beta=b)\n");
      exit(-1);
    }
}

int main(int argc, char ** argv)
{
  double alpha=0.6;
  double beta = alpha;
  //double eps = 1e-8;
  complex double z = 20.0 + 0.0 * I;
  help(argc,argv);
  switch (argc)
    {
    case 2:
      if (argv[1][0]=='t')
	{
	  long i;
	  for(i=0;i<1000000;i++)
	    {
	      mittag_leffler(alpha, beta, z);
	    }
	  return 0;
	}
      if (argv[1][0]=='g')
	{
	  run_grid(alpha,beta);
	  return 0;
	}
      else
      z = atof(argv[1]) + 0.0 * I;
      break;
    case 3:
      z = atof(argv[1]) + 0.0 * I ;
      alpha = atof(argv[2]);
      beta = alpha;
      break;
    case 4:
      z = atof(argv[1]) + 0.0 * I ;
      alpha = atof(argv[2]);
      beta = atof(argv[3]);
      break;
    default:
      break;
    }
  complex double b = mittag_leffler(alpha, beta, z);
  printf("log(mlf(z=%f,alpha=%f,beta=%f)) =  %.8g %+.8gi\n", creal(z), alpha, beta, creal(b), cimag(b));
  return 0;
}
  
#endif /*end standalonemittagleffler*/




