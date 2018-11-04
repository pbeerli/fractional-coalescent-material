#ifndef __SPLINEH__
#define __SPLINEH__
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S P L I N E   R O U T I N E S 
 
 interface part to adaptive spline routines
 from the AMS library
 
 Peter Beerli 1999, Seattle
 beerli@fsu.edu
 
Copyright 1996-2002 Peter Beerli and Joseph Felsenstein, Seattle WA
Copyright 2003-2004 Peter Beerli, Tallahassee FL
 
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

 
$Id: spline.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

extern int dbvssc_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k,
                        long *opt, MYREAL *d0, MYREAL *dnp, MYREAL *d20,
                        MYREAL *d2np, long *constr, MYREAL *eps,
                        MYREAL (*beta) (MYREAL *), MYREAL (*betai) (MYREAL *), MYREAL (*rho) (MYREAL *),
                        MYREAL (*rhoi) (MYREAL *), long *kmax, long *maxstp, long *errc,
                        MYREAL *d, MYREAL *d2, long *diagn, MYREAL *work,
                        long *nwork);

extern int dbvsse_ (MYREAL *x, MYREAL *y, long *np, long *n, long *k,
                        MYREAL *xtab, long *ntab, long *sbopt, long *y0opt,
                        long *y1opt, long *y2opt, long *errc, MYREAL *d,
                        MYREAL *d2, MYREAL *y0tab, MYREAL *y1tab, MYREAL *y2tab,
                        long *erre, MYREAL *work, long *nwork);



#endif
