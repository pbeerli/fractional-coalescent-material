/*------------------------------------------------------
Maximum likelihood estimation
of migration rate  and effectice population size
using a Metropolis-Hastings Monte Carlo [or IS]  algorithm
-------------------------------------------------------


importance sampling for SNP data


Code: Peter Beerli
Model: Michal Palczewski and Peter Beerli

Copyright 2008 Peter Beerli, Tallahassee

    This software is distributed free of charge for non-commercial use
    and is copyrighted. Of course, we do not guarantee that the software
    works and are not responsible for any damage you may cause or have.

*/
#include "migration.h"

void is(world_fmt *world);
void choose_allele(world_fmt *world);

void copy_snp_tips(world_fmt *world, node ** nodep)
{

}

MYREAL
eventtime_snp (world_fmt * world, long thispop, long *ks, char *event, long *whichks)
{
    MYREAL  interval;
    MYREAL  lines;
    MYREAL  denom;
    MYREAL  invdenom;
    MYREAL  rate = world->options->mu_rates[proposal->world->locus];
    MYREAL  invrate = 1./ rate;
    MYREAL  mm = 0.0;
    MYREAL  cc = 0.0;
    long allk = 0;
    for(pop=0; pop < world->numpop; pop++)
      {
	allk += ks[pop];
      }
    
	mm += world->mig0list[pop] * invrate ;
	cc += (ks[pop] * (ks[pop]-1.) * (1.0 / (proposal->param0[pop]*rate)));
      }
    denom    = mm + cc;
    invdenom = 1.0 / denom;
    interval = (-(LOG (UNIF_RANDUM ())) * invdenom) ;
    if ((UNIF_RANDUM ()) < (mm * invdenom))
      {
	*event = 'm';
	return interval;
      }
    else
      {
	*event = 'c';
	return interval;
      }
    error("eventtime_snp failed");
    return -1;
}

void is(world_fmt *world)
{
  node *oldroot;
  node *newroot;
  long *ks;
  long ind;
  node **nodep;
  loong nnodep=0;
  nodep = (node **) mymalloc(world->sumtips,sizeof(node *));
  // create tree
  copy_tree(oldroot,newroot);
  // generate list for easy access
  nnodep = generate_nodelist(newroot,nodep);
  // how many lineages are in each allele class
  ks = (long *) mycalloc(world->numpop * 3, sizeof(long));
  sumks = generate_k(nodep, nnodep, ks);
  // generate a time and event
  while(sumks>0)
    {
      ind = RANDINT(0,sumks); 
      pop = nodep[ind].topop;
      topop = pop;
      age += eventtime_snp(world, &pop, ks, &event);

      if(event == 'c')
	{
	  pick_other(ks,world->numpop, pop, nodep, nnodep);
	  sumks--;	  
	}
      else
	{
	  thenode = nodep[ind];
	  frompop = pick_population(ks, world,thenode);
	  thisnode = add_migration(world, thenode,
				  frompop,
				  topop,
				  age - thenode->tyme);
	}
    }
  calc_like();
}


///
/// generate list of allele counts per populations and total
/// ks is total,allele1,allele2 per population

      
// pick alleletype according to ratio
// if(RANDUM() < binomial(ka,2)/(binomial(ka,2)+binomial(kb2))
// pick two nodes of type x: if in different populations continue
// combine two nodes if coalescence else do nothing but increase time
// calc likelihood  
void pick_other(long *ks, long numpop, long pop, nde *thenode, node **nodep, long nnodep)
{
  // we need to pick a second line with the same population and allele
  // we do this in a two step procedure
  // picking allele
  long j=pop*3;
  long nallelepops;
  double r = RANDUM();
  long rk = ks[j + 1]/ks[j];
  if (r > rk)
    {
      nallele = ks[j+2];
      nallelepops = j+2;
    {
  else
    {
      nallele = ks[j+1];
      nallelepops = j+2;
    }
  // individual
  if (nallele > 1)
    {
      i =RANDOM_INT(0,nallele-1);
      coalesce(thenode, nodep[pop * word->numpop + i]);
    }
}


