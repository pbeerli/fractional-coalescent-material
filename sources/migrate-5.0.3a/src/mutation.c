//
// new mutation model
//
//
/*
 (x) Peter Beerli 2013 Tallahassee FL
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
//data reading
//read_multiple_infile() //put them into the data->yy array (locus x pop x ind)
//adjust_linkage() // allow to move the different loci in/out of linkage groups
//set_mutation_model() // set the model for data type and allow for two different
//                     // types for sequence data and msat data.
//make_values()  // put yy into the tree acknowledging the linkage
//set_inheritance_scalar() // specify the inheritance of each locus
//
//cond_likelihood() // calculate the conditional likelihood over linkage group
//pseudo_cond_like() // same as above but uses fake tree parts
//
#include "migration.h"

// data structures
typedef struct xarray
{
  double **sites;
  long *model;
  long *locuslength;
} xarray;

//data_fmt related
//long **linkage; // for example [(1,3,4,5,7),2,6,(8,9)]
//long *model; // (0,0,0,0,0,1,2,1,2)

//mutation models
//
// allozyme model
//    - infinite allele model
//    - k-allele model
//
//   Q = {{-1 1/n 1/n 1/n ...[n]},{1/n,-1,1/n,1/n,1/n,..[n]}, ..[n]} mu
//   mu = -1/Sum[1/n q_ii,{i,1,n}]
//
// msat model
//    - stepwise mutation model
//    - see Watkins for extensions of the stepwise mutation model
//    - brownian motion model
//
//   Q = {{-Sum[i,{i,1,n}] 1 1/2 1/3 1/4 1/5 ...[n]},{1 -sum[i,{i,1,n}] 1 2 3 4 ...},{2 1 -sum[i,{i,1,n}] 1 2 3 ... }. ...} mu
//   mu = -1 / Sum[1/n q_ii,{i,1,n}]
//
// sequence model
//    - implement the GTR with GAP model and allow for Huelsenbeckian model approach (see Huelsenbeck, Larget, Alfaro)
//
//
// Q = {{a11,         rac pic,   rag pig,   rat pit,   ragap pigap},
//      {rac pia,     a22,       rcg pig,   rct pit,   rcgap pigap},
//      {rag pia,     rcg pic,   a33,       rgt pit,   rggap pigap},
//      {rat pia,     rct pic,   rgt pig,   a44,       rtgap pigap},
//      {ragap pia,   rcgap pic, rggap pig, rtgap pit, a55}}
//
// a11 = -  rac pic - rag pig - rat pit - ragap pigap
// a22 = - rac pia - rcg pig - rct pit - rcgap pigap
// a33 = -rag pia - rcg pic - rgt pit - rggap pigap
// a44 = - rat pia - rct pic - rgt pig - rtgap pigap
// a55 = - ragap pia - rcgap pic - rggap pig - rtgap pit
//
// mu = -1/(pia a11 + pic a22 + pig a33 + pit a44 + pigap a55)
//
// T = MatrixExp[t Q mu] with t the branchlength, this is of course rate*t' to accomodate the different rate categories
//
// model = expressed as a sequence of values that set the rij values equal,
// model = rac rag rat rcg rct rgt raGap rcGap rgGap rtGap
//         for example model= 1111110000 is equal to F81
//         model=1234560000 is the GTR model
//         model=123456789a is the GTR + Gap model
// 
// datastorage device: base frequencies , rate parameters
// represent the known models with additional calculation for the the printout but
// this is not needed for the calculations (I believe if we have the transition probs
// the likelihood is calculated as
// per site
// condL[gi_, gj_, ti_, tj_] := Module[{},
//    pi = p[ti];
//    pj = p[tj];
//    hi =  Plus @@ (pi gi);
//    hj = Plus @@ (pj gj);
//    hi hj 
//    ]
// for the whole sequence


//  the vectors qmatrix and rmatrix is in fact a {{x1,x2,x3,...},{...}, ...} matrix
void MatrixExp(double *rmatrix, double *qmatrix, double qsize, double t)
{
  long i;
  const long qsize2 = qsize * qsize; //this assumes that the qmatrix and rmatrix are contigously aligned
  double zfact = 1;
  double z = 0.;
  double tz;
  double sum;
  double oldsum = HUGE;
#ifdef ROUGH
  const double compare = 0.001;
#else
  const double compare = EPSILON;
#endif
  double err = HUGE;
  while(err<compare && z < 1000)
    {
      z += 1.0;
      zfact *= z;
      sum = 0.0;
      tz = t / zfact;
      for(i=0;i<qsize2; i++)
	{
	  *(qmatrix+i) += (*(rmatrix+i)) * tz;
	  sum += *(qmatrix+i);
	}
      err = fabs(sum - oldsum);
      oldsum = sum;
    }
}

void  Both_Qmatrix_to_Tmatrix(double t1, double t2, double **q1, double **q2, double **Rmatrix, long Qsize, long number_of_qmatrices)
{
  long i;
  for(i=0;i<number_of_qmatrices;i++)
    {
      MatrixExp(Rmatrix[i],q1[i],Qsize, t1);
      MatrixExp(Rmatrix[i],q2[i],Qsize, t2);
    }
}

///
/// conditional likelihood calculator
/// takes times conditional likelihoods, a vector of  Q-matrix and its individual size, 
/// a vector of models, the number of models and the 
double condL(double t1, double t2, double *cLike1, double *cLike2, double *cLike3, 
	     double **Rmatrix, long Qsize, long number_of_qmatrices, long *model, long modelLength)
{
  long size;
  long nuc;
  long snuc;
  long pos;
  long site;
  long nQ;
  // fill these with the transition probabilities
  double *q1 = get_matrix(Qsize, modelLength);
  double *q2 = get_matrix(Qsize, modelLength);
  double *tr1;
  double *tr2;
  double h1;
  double h2;

  Both_Qmatrix_to_Tmatrix(t1, t2, q1, q2, Rmatrix, Qsize, number_of_qmatrices);
  
  for(site=0; site<modelLength; site++)
    {

      tr1 = q1 + model[site];
      tr2 = q2 + model[site];

      for(nuc=0;nuc<Qsize;nuc++)
	{
	  h1 = 0.0;
	  h2 = 0.0;
	  snuc = site + nuc;
	  nQ = nuc*Qsize;
	  for(pos=0;pos<Qsize;pos++)
	    {
	      h1 += cLike1[snuc] * t1[nQ + pos] ;
	      h2 += cLike2[snuc] * t2[nQ + pos] ;
	    } 
	  cLike3[snuc] = h1 * h2;
	}
    }
}


void calculate_condlike(world_fmt *world, double **cLike1, double **cLike2, double **cLike3)
{
  long r;
  double rate;
  mu_fmt *mumodel = world->mumodel;
  for(r=0;r<mumodel->rcategs;r++)
    {
      rate = mumodel->rates[r];
      condL(rate * t1, rate * t2, cLike1, cLike2, cLike3, 
	    mumodel->Qmatrix, mumodel->QSize, mumodel->number_of_qmatrices, 
	    mumodel->seqmodel, mumodel->seqmodelLength)
    }
}

double treelikelihood(world_fmt *world)
{
  double loglike;
  long a;
  MYREAL term = 0.0;
  node *nn = crawlback (world->root->next);
  set_dirty (nn);
  smooth (world->root->next, crawlback(world->root->next), world, world->locus);
  adjustroot (world->root);
  for(site=0; site < mumodel->modelLength; site++)
    {
      calc_rate_site(cLike1,basefrequencies, tjs);
      ts = calc_mean_rates_site(tjs, rateprobs);
      logts = log(ts) + scaler[s];
      calc_relative_rates_site(tjs,ts, cjs);
      sumS += ws * logts; // what is ws (text says w_s weighted sume of rate averages
    }
  memset(like,0,sizeof(double)*mumodel->likelength);
  // order of calculation needs to follow the real sequence 
  // lambda is the patch-length rate size=2-> lambda=1/2
  for(site=0; site < mumodel->totalDataLength; site++)
    {
      sumScs = lambda * calc_mean_rates_site(like,rateprobs);
      recalc_like(cjs, lambda, sumScs, like , newlike);
      swap(like,newlike);
    }
  loglike = sumS + log(calc_mean_rates_site(like,rateprobs))
  return loglike;
}

void nuview_gap(double *mother, double *left, double *right, scaler)
{
  long u;
  long i;
  long q;
  double *q1 = get_matrix(Qsize, modelLength);
  double *q2 = get_matrix(Qsize, modelLength);
  double *tr1;
  double *tr2;
  double h1;
  double h2;

  Both_Qmatrix_to_Tmatrix(t1, t2, q1, q2, Rmatrix, Qsize, number_of_qmatrices);

  for(u=0;u<endsite;u++)
    {
      double *ll = left->seq[u];
      double *rr = right->seq[u];
      double *mm = mother->seq[u];
      
      for(i=0;i<K+1;i++)
	{
	  mm[i] = 0.0;
	  for(q=0;q<K;q++)
	    {
	      mm[i] += (ll[q] * trans_prob[i][q] + (is_gap(ll[q]) * trans_prob[i][GAP]))
		*(rr[q] * trans_prob[i][q] + (is_gap(right[q]) * trans_prob[i][GAP]));
	    }
	}
    }
}

// use this in smooth for the to calculate the extra normalizing material
// leaves are already initialized for star=1 this is the same for all sites
MYINLINE void calc_prob_star(double left_star, double right_star, double *mother_star, 
		    double x1m__t_left, double x1m__t_right)
{
  *mother = left_star * x1m_t_left * right_star * x1m_t_right; 
}


// calculates 1 - xi_t
// 1 - lambda/sum * (1.0 - exp(sum * t));
MYINLINE double xi1m_t(double lambda_sumratio, double lambdasum, double t)
{
  return 1.0 - lambda_sumratio + (1.0 - exp(lambdasum * t));
}



// calculate l1pop_gap
// using scaler to make sure that the numbers do not get too small
// once in logs we are safe, therefore turn unscaled log value;
MYINLINE double calc_l1prob_gap_TRp(double prob_gap_root_scaled, double scaler)
{
  return log(1.0 - prob_gap_root_scaled * exp(scaler)
}

// calculate lprob_star/*TRp*/
// using scaler to make sure that the numbers do not go too small
// once in logs we are safe, therefore turn unscaled log value;
double calc_prob_star_TRp(double prob_star_root_scaled, double scaler, double p)
{
  return log( (1.0 - p) * prob_star_root_scaled) + scaler;
}

// calculate the tree likelihood assuming all sites are indepdendent, 
// and no site rates;
// l1prob_gap can be calcualted using 
// log1p(-prob_gap)
MYINLINE double calc_treelike(double lprob_star/*_TRp*/, double l1prob_gap/*_TRp*/, 
		     double )
{
  double lnL=0.;
  lnL = lprob_star - l1prob_gap;
  for(u = 0; u < endsites; u++)
    {
      lnL += rootLnL[u];
    }
  returm lnL;
}
