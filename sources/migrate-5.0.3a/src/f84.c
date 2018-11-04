/*------------------------------------------------------
 Bayesian Inference and Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 F84 withOUT GAPS
 
 
 Peter Beerli 2008, Tallahassee
 beerli@fsu.edu
 
 Copyright 2008-present Peter Beerli
 
 This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
 $Id:$
 
-------------------------------------------------------*/
/* \file f84.c

*/
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "watterson.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif

/* prototypes ------------------------------------------- */

///
/// condlike will replace nuview
void condlike(world_fmt *world)
{

  seq_fmt * seq = world->seq;

  const long endsite = seq->endsite;
  register float lw1 = -b1 * seq->fracchange;
  register float lw2 = -b2 * seq->fracchange;
  register float freq[8];//a,c,g,t,ar,gr,cy,ty

  freq[0] = seq->freqa;
  freq[1] = seq->freqc;
  freq[2] = seq->freqg;
  freq[3] = seq->freqt;
  freq[4] = seq->freqar;
  freq[5] = seq->freqgr;
  freq[6] = seq->freqcy;
  freq[7] = seq->freqty;
  for (j=0; j < mutmodels; j++)
    {
      (*setup_calculator)[j](world,&freq, b1, b2);
    }
  for (i = 0; i < endsite; i++)
    {
      (*site_calculator)[i](xx1,xx2,xx3, world, freq, i);
    }  
}

static void * (*site_calculator) (node *, node *, node *, world_fmt *, long );
///
/// sets up the conditional likelihood calculator based on the mutation model used
/// the conditional likelihood calculator is mapped for each "site"
setup_site_calculator(long endsite)
{
  long i;
  if(site_calculator == NULL)
    {
      site_calculator = (void *) malloc(sizeof(void) * endsite);
    } 
  else
    {
      site_calculator = (void *) realloc(site_calculator, sizeof(void) * endsite);
    }
  for (i = 0; i < endsite; i++)
    {
      switch(datatype[i])
	{
	case F84:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) f84_calculator;
	  break;
	case F84GAP:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) f84gap_calculator;
	  break;
	case BROWN:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) brownian_calculator;
	  break;
	case MSAT:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) msat_calculator;
	  break;
	case ALLOZYME:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) allozyme_calculator;
	  break;
	case TWOSTATE:
	  site_calculator[i] = (void (*)(node *, node *, node *, world_fmt *, long )) twostate_calculator;
	  break;
	}
    }
}

void f84gap_calculator(node * xx1,node * xx2, node * xx3, world_fmt * world, long datatype)
{
  const struct mutationmodel = world->seq->mutationmodel[datatype];

}




setup_sites()
{
  // fill loci into containers
  // 
}

// calculator
yy1 = 1 - EXP (tbl00->ratxv * lw1);
yy2 = 1 - EXP (tbl00->ratxv * lw2);

T1 = freq * xx1; // each element
T2 = freq * xx2;

sum(T1)*yy1
sum(T2)*yy2


 
