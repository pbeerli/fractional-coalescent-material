#ifndef COMBRO_INCLUDE
#define COMBRO_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice world size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 C O M B I N E   L O C I   R O U T I N E S 
 
 
 
 Peter Beerli 1996, Seattle
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



$Id: combroyden.h 2157 2013-04-16 22:27:40Z beerli $
 
-------------------------------------------------------*/

#include "migration.h"

extern long estimateParameter (long rep, long G, world_fmt * world,
                                   option_fmt * options, MYREAL **dd, long chain,
                                   char type, long kind, long repkind);

extern MYREAL calc_loci_like (helper_fmt * helper, MYREAL *param,
                                  MYREAL *lparam);

extern void calc_apgg (MYREAL thetax, MYREAL **apgg, nr_fmt * nr,
                           timearchive_fmt * atl, long which);
extern void set_gamma_param (MYREAL *paramn, MYREAL *paramo, MYREAL *lparamn,
                                 MYREAL *lparamo, MYREAL theta, nr_fmt * nr);
extern void calc_apg0 (MYREAL *apg0, nr_fmt * nr, timearchive_fmt * tyme,
                           MYREAL thetax, MYREAL theta1);
extern void which_calc_like (long repkind);
extern void calc_loci_param (nr_fmt * nr, MYREAL *lparam, MYREAL *olparam,
                                 MYREAL *dv, MYREAL lamda, long nnn);
extern void set_replicates (world_fmt * world, long repkind, long rep,
                                long *repstart, long *repstop);
extern void prepare_broyden (long kind, world_fmt * world,
                                 boolean * multilocus);
extern void simple_loci_derivatives (MYREAL *d, nr_fmt * nr,
                                         timearchive_fmt ** tyme, long locus);
extern void copy_and_clear_d (nr_fmt * nr);
extern void add_back_d (nr_fmt * nr);
extern void setup_parameter0_standard (world_fmt * world, nr_fmt * nr,
                                           long repkind, long repstart,
                                           long repstop, long loci, long kind,
                                           boolean multilocus);

extern void combine_gradient (nr_fmt * nr, helper_fmt * helper, MYREAL *gxv);

extern long maximize (MYREAL **thisparam, world_fmt * world, 
		      nr_fmt * nr,
		      MYREAL **hess,
                          long analystype, long repkind);
extern void set_logparam (MYREAL *loga, MYREAL *a, long size);
extern void set_expparam (MYREAL *expa, MYREAL *a, long size);
extern long do_profiles (world_fmt * world, nr_fmt * nr, MYREAL *likes,
                             MYREAL *normd, long kind, long rep, long repkind);
extern void fill_helper (helper_fmt * helper, MYREAL *param, MYREAL *lparam,
                             world_fmt * world, nr_fmt * nr);

extern void replace_with (int mode, MYREAL *a, MYREAL *b, MYREAL *c, MYREAL m,
                              MYREAL *la, MYREAL *lb, MYREAL *lc, MYREAL ll);
MYREAL quadratic_lamda (MYREAL lamda, MYREAL a, MYREAL b, MYREAL c, MYREAL la,
                        MYREAL lb, MYREAL lc);


#endif /* COMBINE_INCLUDE */
