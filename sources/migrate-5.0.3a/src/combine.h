#ifndef COMBINE_INCLUDE
#define COMBINE_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice world size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 C O M B I N E L O C I (NEWTON RAPHSON)  R O U T I N E S 

 combines loci by estimating a gamma shape parameter
 

 Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
 (c) Peter Beerli Tallahassee 2013
 
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
 

 $Id: combine.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

#include "migration.h"

extern long combine_loci (world_fmt * world);
//extern void create_nr(nr_fmt *nr,long numpop, long G);
//extern void destroy_nr(nr_fmt *nr);
extern double calc_loci_like (nr_fmt * nr, timearchive_fmt * atl, long loci,
                                  boolean boolgamma);
extern void calc_gamma (nr_fmt * nr);
extern void calc_apgg (double thetax, double **apgg, nr_fmt * nr,
                           timearchive_fmt * atl, long which);
extern void set_gamma_param (double *paramn, double *paramo, double theta,
                                 nr_fmt * nr);
extern long broyden_driver (timearchive_fmt * tyme, long loci,
                                world_fmt * world, double **covariance,
                                char ***plane, long *profiles, double *value,
                                long profilenum, double *profparam,
                                double *proflike, double *profnorm);
extern void calc_apg0 (double *apg0, nr_fmt * nr, timearchive_fmt * tyme,
                           double thetax, double theta1);
#endif /* COMBINE_INCLUDE */
