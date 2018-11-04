/* Cholesky update

   updating  L . D . L' 

   using Gil et al.


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




   $Id: cholesky.c 2067 2012-07-27 20:59:32Z beerli $
 */
#include "migration.h"
#include "cholesky.h"


/* functions profiles */
void update_B (double **l, double *d, double lamda, double *dv, double *gamma,
               double *delta, long n, double **work);
void update_cholesky (double **l, double *d, double *v, double sign,
                      double **vv, double *beta, double *t, double *p,
                      long n);
void calc_v_w (double *v, double *w, double *sign1, double *sign2,
               double lamda, double *dv, double *gamma, double *delta,
               long n);
void calc_direction (double **l, double *d, double *dv, double *gxv,
                     double **ll, long n);
void calc_lower_transpose (double **l, double *d, double **ll, long n);
/*double brent (double ax, double bx, double cx, double (*func) (),
   double tol, double *xmin, void *helper); */
void solve_p (double *p, double **l, double *d, double *v, double *psquare,
              long n);
/*
   B = B0 + vv'
   B = LDL' + vv'
   = L(D+pp')L'
   where
   p is the solution to Lp=v
   . .
   B = LXL'    X= (LDL')
   = LDL'

 */

/* update the cholesky factors
   Gil, P.E., W. Murray, and M.H. Wright. 1981
   Practical optimization, Academic Press

   updating the factors of B = L D L

   L     n x n lower diagonal matrix
   D     n x 1 diagonal matrix 
   LAMBDA    1 length to jump into search direction
   DV    n x 1 search direction
   GAMMA n x 1 differences of first derivatives g_i - g_(i-1)
   DELTA n x 1 differences of paramter values x_i - x_(i-1)
   N         1 number of elements n
   WORK  (2n+4) x n work array

   Peter Beerli 1998 (beerli@fsu.edu)

 */

void
update_B (double **l, double *d, double lamda, double *dv, double *gamma,
          double *delta, long n, double **work)
{
    double signv, signw;
    double *v = work[0];
    double *w = work[1];
    double *beta = work[2];
    double *t = work[3];
    double *p = work[4];
    double **vv = work + 5;
    calc_v_w (v, w, &signv, &signw, lamda, dv, gamma, delta, n);
    update_cholesky (l, d, v, signv, vv, beta, t, p, n);
    update_cholesky (l, d, w, signw, vv, beta, t, p, n);
}


/* are the shortcuts correct, compare to paper */
void
update_cholesky (double **l, double *d, double *v, double sign, double **vv,
                 double *beta, double *t, double *p, long n)
{
#ifdef MYREAL == float
  const MYREAL eps = FLT_EPSILON ;
#else
  const  MYREAL eps = DBL_EPSILON ;
#endif
    long j, r;
    double psquare = 0;
    if (sign > 0.0)
    {
        t[0] = 1.;
        memcpy (vv[0], v, sizeof (double) * n);
        for (j = 0; j < n; j++)
        {
            p[j] = v[j];
            t[j + 1] = t[j] + p[j] * p[j] / d[j];
            beta[j] = p[j] / (d[j] * t[j + 1]);
            d[j] *= t[j + 1] / t[j];
            for (r = j + 1; r < n; r++)
            {
                vv[j + 1][r] = vv[j][r] - p[j] * l[r][j];
                l[r][j] += beta[j] * vv[j + 1][r];
            }
        }
    }
    else
    {
        solve_p (p, l, d, v, &psquare, n);
        t[n] = 1. - psquare;
        if (t[n] <= eps)
            t[n] = eps;
        for (j = n - 1; j >= 0; j--)
        {
            t[j] = t[j + 1] + p[j] * p[j] / d[j];
            beta[j] = -p[j] / (d[j] * t[j + 1]);
            vv[j][j] = p[j];
            d[j] *= t[j + 1] / t[j];
            for (r = j + 1; r < n; r++)
            {
                vv[j][r] = vv[j + 1][r] + p[j] * l[r][j];
                l[r][j] = l[r][j] + beta[j] * vv[j + 1][r];
            }
        }
    }
}

void
solve_p (double *p, double **l, double *d, double *v, double *psquare, long n)
{
    long i, r;
    p[0] = l[0][0] / v[0];
    for (i = 0; i < n; i++)
    {
        p[i] = v[i];
        for (r = 0; r < i; r++)
        {
            p[i] += -l[r][i] * p[r];
        }
        p[i] /= l[i][i];
        *psquare += p[i] * p[i] / d[i];
    }
}

void
calc_v_w (double *v, double *w, double *sign1, double *sign2, double lamda,
          double *dv, double *gamma, double *delta, long n)
{
    long i;
    double denom1 = 0.;
    double denom2 = 0.;

    for (i = 0; i < n; i++)
    {
        denom1 += gamma[i] * dv[i];
        denom2 += lamda * delta[i] * gamma[i];
    }
    *sign1 = denom1 < 0 ? -1. : 1;
    *sign2 = denom2 < 0 ? -1. : 1;
    for (i = 0; i < n; i++)
    {
        v[i] = gamma[i] / sqrt (fabs (denom1));
        w[i] = delta[i] / sqrt (fabs (denom2));
    }
}


/* solve LDL' p = -g */

void
calc_direction (double **l, double *d, double *dv, double *gxv, double **ll,
                long n)
{
    long i, j;
    double sum = 0.0;
    long nn = n;

    calc_lower_transpose (l, d, ll, nn);

    for (i = 0; i < nn; i++)
    {
        sum = gxv[i];
        for (j = 0; j < i; j++)
        {
            sum -= ll[i][j] * (gxv[j]);
        }
        dv[i] = sum;
    }
    for (i = nn - 1; i >= 0; i--)
    {
        sum = dv[i];
        for (j = i + 1; j < nn; j++)
        {
            sum -= ll[i][j] * dv[j];
            dv[i] = sum / ll[i][i];
        }
    }
    for (i = 0; i < nn; i++)
    {
        if (fabs (dv[i]) > 100000.)
            dv[i] = (dv[i] < 0 ? -100000. : 100000);
    }
}


void
calc_lower_transpose (double **l, double *d, double **ll, long n)
{
    long i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            ll[i][j] = l[i][j] * sqrt (d[i]);
            ll[j][i] = ll[i][j];
        }
        ll[i][i] = l[i][i] * sqrt (d[i]);
    }
}


/*
   #define ITMAX 100
   #define CGOLD 0.3819660
   #define ZEPS 1.0e-10

   #define SIGN(a,b)  ((b) > 0.0 ? fabs(a) : -fabs(a))
   #define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

   double 
   brent (double ax, double bx, double cx, double (*func) (),
   double tol, double *xmin, void *helper)
   {
   long iter;
   double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w,
   x, xm;
   double e = 0.0;

   a = ((ax < cx) ? ax : cx);
   b = ((ax > cx) ? ax : cx);
   x = w = v = bx;
   fw = fv = fx = (*func) (x, helper);
   for (iter = 1; iter <= ITMAX; iter++)
   {
   xm = 0.5 * (a + b);
   tol2 = 2.0 * (tol1 = tol * fabs (x) + ZEPS);
   if (fabs (x - xm) <= (tol2 - 0.5 * (b - a)))
   {
   *xmin = x;
   return fx;
   }
   if (fabs (e) > tol1)
   {
   r = (x - w) * (fx - fv);
   q = (x - v) * (fx - fw);
   p = (x - v) * q - (x - w) * r;
   if (q > 0.0)
   p = -p;
   q = fabs (q);
   etemp = e;
   e = d;
   if (fabs (p) >= fabs (0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
   d = CGOLD * (e = (x >= xm ? a - x : b - x));
   else
   {
   d = p / q;
   u = x + d;
   if (u - a < tol2 || b - u < tol2)
   d = SIGN (tol1, xm - x);
   }
   }
   else
   {
   d = CGOLD * (e = (x >= xm ? a - x : b - x));
   }
   u = (fabs (d) >= tol1 ? x + d : x + SIGN (tol1, d));
   fu = (*func) (u, helper);
   if (fu <= fx)
   {
   if (u >= fx)
   a = u;
   else
   b = u;
   SHFT (v, w, x, u);
   SHFT (fv, fw, fx, fu);
   }
   else
   {
   if (fu <= fw || w == x)
   {
   v = w;
   w = u;
   fv = fw;
   fw = fu;
   }
   else if (fu <= fv || v == x || v == w)
   {
   v = u;
   fv = fu;
   }
   }
   }
   FPRINTF (stderr, "too many iterations in brent\n");
   *xmin = x;
   return fx;
   }

 */
