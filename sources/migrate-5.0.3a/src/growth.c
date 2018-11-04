 /*------------------------------------------------------
 inference of population parameters
 using a Metropolis-Hastings Monte Carlo algorithm
 -------------------------------------------------------
 growth routines
 
 Peter Beerli 2013, Tallahassee
 beerli@fsu.edu
 
 Copyright 2017 Peter Beerli
 
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
 
*/
#include "growth.h"
#include "sighandler.h"
#include "bayes.h"
void reset_growth(world_fmt * world);
#if defined(MPI) && !defined(PARALIO) /* */
void print_growth_record(float *temp, long *z, world_fmt *world);
#else
void  print_growth_record(char *temp, long *c, world_fmt * world);
#endif
boolean init_growpop(worldoption_fmt * wopt, option_fmt *options, long numpop);


void reset_growth(world_fmt * world)
{
  long i;
  if (world->has_growth)
    {
      for(i=0;i<world->options->growpops_numalloc;i++)
	{
	  world->growth[world->options->growpops[i]-1]=0.0;
	}
    }
}


boolean init_growpop(worldoption_fmt * wopt, option_fmt *options, long numpop)
{
  long i;
  boolean use_growth = FALSE;
  for (i=0;i<options->growpops_numalloc;i++)
    {
      if (options->growpops[i]>0)
	{
	  use_growth = TRUE;
	  break;
	}
    }
  if (use_growth)
    {
      wopt->growpops = (long*) mycalloc(numpop,sizeof(long));
      memcpy(wopt->growpops,options->growpops, sizeof(double) * (size_t) options->growpops_numalloc);
      wopt->growpops_numalloc = options->growpops_numalloc;
    }
  else
    {
      wopt->growpops = NULL;
      wopt->growpops_numalloc = 0;
    }
  return use_growth;
}

void init_growth(world_fmt * world, long numpop)
{
  long i;
  world->grownum = 0;
  if(world->options->growpops != NULL)
    {
      world->has_growth = TRUE;
      world->growth = (double*) mycalloc(numpop,sizeof(double));
      world->savegrowth = (double*) mycalloc(numpop,sizeof(double));
      for (i=0; i < world->options->growpops_numalloc; i++)
	{
	  long x = world->options->growpops[i];
	  if (x != 0)	   
	    world->growth[x-1] = 1.0;
	}
      world->grownum = 0;
      for (i=0; i < numpop; i++)
	{
	  if(world->growth[i] > 0.0)
	    {
	      world->grownum += 1;
	    }
	  world->growth[i]=0.0;
	}
    }
  else
    {
      world->has_growth = FALSE;
    }
}


#if defined(MPI) && !defined(PARALIO) /* */
void print_growth_record(float *temp, long *z, world_fmt *world)
{
  long i;
  if (world->has_growth)
    {
      for (i=0; i<world->numpop;i++)
	{
	  temp[(*z)++] = world->growth[i];
	}
    }
}
#else /*not MPI or MPI & PARALIO*/
void  print_growth_record(char *temp, long *c, world_fmt * world)
{
  long i;
  if (world->has_growth)
    {
      for (i=0; i<world->numpop;i++)
	{
	  *c += sprintf(temp+ *c,"\t%f", world->growth[i]); 
	} 
    }
}
#endif

void construct_locusgrowth_histogram(world_fmt *world, long locus, MYREAL *mini, MYREAL *maxi, double **results)
{
  bayes_fmt *bayes = world->bayes;
  long i;
  long j0;
  double themean;
  double thestd;
 
  long numbin = 0;
  long pa;
 
  long rpa;
  long total=0;
  long *bins = bayes->histogram[locus].bins;
  boolean *visited;
  
  long grownum = world->options->growpops_numalloc;
  long np = world->numpop2 + bayes->mu + 2 * world->species_model_size;
  
  visited = (boolean *) mycalloc(world->options->growpops_numalloc, sizeof(boolean));
  for(i=0;i<np;i++)
    numbin += bins[i];
  
  for(j0=0; j0 < grownum; j0++)
    {
      long pick = world->options->growpops[j0];
      if (pick == 0)
	continue;
      else
	rpa = pick-1;
      if(!visited[rpa])
        {
	  themean = 0.0;
	  thestd = 0.0;
	  total = 0;
	  pa = rpa + np; //this points into the histograms (!! all other parameters, and then growth!!)	  
	  //np=#all param before growth, rpa= #growhindex, numbin=#bins before #grwothindex
	  construct_param_hist(world,locus,np, rpa,numbin, mini, maxi, results, &total,&themean,&thestd);
	  world->bayes->histogram[locus].means[pa] = themean/total;
	  world->bayes->histogram[locus].stds[pa]  = thestd / total;
	  world->bayes->histtotal[locus*np+pa] = (MYREAL) total;
	  visited[rpa] = TRUE;
	  numbin += bins[pa];
        }
    }
}

  
