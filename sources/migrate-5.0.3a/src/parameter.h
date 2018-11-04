





#ifndef PARAM_INCLUDE
#define PARAM_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice world size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 P A R A M E T E R E S T I M A T I O N   R O U T I N E S 

 estimates parameter for each locus
 using a Newton-Rapshon maximization
 

 (c) Peter Beerli 1996, Seattle
 beerli@fsu.edu
 
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

 $Id: parameter.h 2157 2013-04-16 22:27:40Z beerli $
-------------------------------------------------------*/

extern void estimateParameter (timearchive_fmt * tyme, long G,
                                   world_fmt * world, double **dd, long chain,
                                   char type, char **plane);
extern double probG (double *param, tarchive_fmt * tl, long numpop);
extern double logprob_noevent (world_fmt * world, long interval);
extern double sum_migprob (world_fmt * world, long pop, long interval);
/* calculates the likelihood for the new parameter set */
extern double calc_like (nr_fmt * nr, tarchive_fmt * tyme, long G);
extern void reset_nr (nr_fmt * nr);
extern void derivatives_to_logderivatives (nr_fmt * nr);
/* calculate the norm sqrt(sum(v*v)) */
extern double norm (double *d, long size);
extern boolean is_singular (double **dd, long n);
extern void calc_param (nr_fmt * nr, double *param, double lamda);
extern void param_all_adjust (nr_fmt * nr, double *param, long gamma_param);
extern void derivatives (long trials, nr_fmt * nr, tarchive_fmt * tl, long G,
                             double *param, boolean forloci);
extern void derivatives_to_logderivatives (nr_fmt * nr);
extern void calc_cov (double **dd, double *d, double *param, long n);
extern void free_nr (nr_fmt * nr);
extern double vector_max (double *v, long size);
#endif /*PARAM_INCLUDE */
