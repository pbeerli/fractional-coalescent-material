/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 S E Q U E N C E S   R O U T I N E S 
 F84 with GAPS
 
 
 Peter Beerli 2[0]7, Tallahassee
 beerli@fsu.edu
 
 Copyright 2[0]7 Peter Beerli
 
  This software is distributed free of charge for non-commercial use
 and is copyrighted. Of course, we do not guarantee that the software
 works and are not responsible for any damage you may cause or have.
 
 
 $Id:$
 
-------------------------------------------------------*/
/* \file f84gap.c

*/
#ifdef GAP
#include "migration.h"
#include "sighandler.h"
#include "tools.h"
#include "migrate_mpi.h"
#include "watterson.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */



MYINLINE MYREAL f1(const MYREAL v, const MYREAL alpha, const MYREAL gamma)
{
  return EXP(v * (gamma-alpha));
}


MYINLINE MYREAL f2(const MYREAL v, const MYREAL gamma)
{
  return EXP(v * (-gamma));
}


// f3 is same function as f2 but call with alpha instead of gamma

MYINLINE MYREAL p3(const MYREAL thisf2, const MYREAL pi_n)
{
  return pi_n * (1.0 - thisf2);
}

MYINLINE MYREAL p4(const MYREAL thisf2, const MYREAL thisp3)
{
  return thisf2 + thisp3;
}

MYINLINE MYREAL p2(const MYREAL pi_n_pi_cl_ratio, const MYREAL thisf2, const MYREAL thisf3, const MYREAL thisp3)
{
  return thisf2 * (1.0 - thisf3) * pi_n_pi_cl_ratio + thisp3; 
}

MYINLINE MYREAL p1(const MYREAL thisf1, const MYREAL thisp2)
{
  return thisf1 + thisp2;
}

MYINLINE MYREAL f84gap(const long i, const long j, const MYREAL *trprob)
{
  return trprob[i * 5 + j];
}

MYINLINE void setup_f84gap(MYREAL *trprob, const MYREAL v,const MYREAL alpha,const MYREAL gamma,const MYREAL pi_A, const MYREAL pi_C, const MYREAL pi_G, const MYREAL pi_T, const MYREAL pi__, const MYREAL inv_pi_R, const MYREAL inv_pi_Y)
{
  
  const MYREAL thisf1  = f1(v, alpha, gamma);
  const MYREAL thisf2  = f2(v, gamma);
  const MYREAL thisf3  = f2(v, alpha);

  const MYREAL thisp3A = p3(thisf2,pi_A); 
  const MYREAL thisp3C = p3(thisf2,pi_C); 
  const MYREAL thisp3G = p3(thisf2,pi_G); 
  const MYREAL thisp3T = p3(thisf2,pi_T); 
  const MYREAL thisp3_ = p3(thisf2,pi__); 

  const MYREAL thisp2A  = p2(pi_A * inv_pi_R, thisf2, thisf3, thisp3A);
  const MYREAL thisp2C  = p2(pi_C * inv_pi_Y, thisf2, thisf3, thisp3A);
  const MYREAL thisp2G  = p2(pi_G * inv_pi_R, thisf2, thisf3, thisp3A);
  const MYREAL thisp2T  = p2(pi_T * inv_pi_Y, thisf2, thisf3, thisp3A);

  const MYREAL thisp4A = p4(thisf2,thisp3A); 
  const MYREAL thisp4C = p4(thisf2,thisp3C); 
  const MYREAL thisp4G = p4(thisf2,thisp3G); 
  const MYREAL thisp4T = p4(thisf2,thisp3T);
  const MYREAL thisp4_ = p4(thisf2,thisp3_);

  const MYREAL thisp1A = p1(thisf1,thisp2A); 
  const MYREAL thisp1C = p1(thisf1,thisp2C); 
  const MYREAL thisp1G = p1(thisf1,thisp2G); 
  const MYREAL thisp1T = p1(thisf1,thisp2T);


  trprob[0] = thisp1A;
  trprob[6] = thisp1C;
  trprob[12] = thisp1G;
  trprob[18] = thisp1T;
  trprob[24] = thisp4_;

  trprob[1] = thisp3A;
  trprob[2] = thisp2A;
  trprob[3] = thisp3A;
  trprob[4] = thisp3A;

  trprob[5] = thisp2C;
  trprob[7] = thisp3C;
  trprob[8] = thisp2C;
  trprob[9] = thisp3C;

  trprob[10] = thisp2G;
  trprob[11] = thisp3G;
  trprob[13] = thisp3G;
  trprob[14] = thisp3G;

  trprob[15] = thisp3T;
  trprob[16] = thisp2T;
  trprob[17] = thisp3T;
  trprob[19] = thisp3T;

  trprob[20] = thisp4A;
  trprob[21] = thisp4C;
  trprob[22] = thisp4G;
  trprob[23] = thisp4T;
}


///
/// accurate version of conditional likleihood calculator using a scaler
/// and the F84+GAP model
void f84gap_calculator(node * xx1,node * xx2, node * xx3, world_fmt * world, long datatype)
{
  const struct mutationmodel = world->seq->mutationmodel[datatype];
}

void nuview_f84gap_slow (node * mother, world_fmt * world, const long locus)
{
    static long count = 0;
    long i, j, k;
    MYREAL transprobr[9 * 9 * 5 * 5];
    MYREAL transprobq[9 * 9 * 5 * 5];
    MYREAL *tblr;
    MYREAL *tblq;    
    MYREAL gamma;
    MYREAL alpha;
    MYREAL vr;
    MYREAL vq;
    MYREAL xx1t0, xx1t1,xx1t2, xx1t3, xx1t4;
    MYREAL xx2t0, xx2t1, xx2t2, xx2t3, xx2t4;
    MYREAL sa1,sa2,sc1,sc2,sg1,sg2,st1,st2,s_1,s_2;

    MYREAL lw1, lw2;
    //, yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vzsumr1,
    //    vzsumr2, vzsumy1, vzsumy2, sum1, sum2, sumr1, sumr2,
    //    sumy1, sumy2, ww1, ww2, zz1, zz2;
    MYREAL *sxx1 = NULL, *sxx2 = NULL;
    MYREAL sxx3m, tempsxx3m;
    node *q, *r;
    sitelike *xx1, *xx2, *xx3;
    long rcategs = world->options->rcategs;
    long categs = world->options->categs;
    tbl_fmt tbl = world->tbl;
    MYREAL freqa, freqc, freqg, freqt, freq_, freqr, freqy;
    seqmodel_fmt *seq;
    //valrec *tbl[0][0];
    valrec *tblij;
    valrec *tbljk;
    seq = world->data->seq;
    freqa = seq->freqa;
    freqc = seq->freqc;
    freqg = seq->freqg;
    freqt = seq->freqt;
    freqr = seq->freqr;
    freqy = seq->freqy;
    freq_ = seq->freq_;
    q = crawlback (mother->next);
    r = crawlback (mother->next->next);
    vq = q->v;
    vr = r->v;
    sxx1 = q->scale;
    sxx2 = r->scale;
    //    tbl[0] = tbl[0];
    alpha = seq->xi;
    gamma = seq->xv;

    if (rcategs == 1 && categs == 1)
    {
      tblr = &transprobr[0];
      tblq = &transprobq[0];
      setup_f84gap(tblq, vq , alpha, gamma, 
		   freqa, freqc, freqg, freqt, freq_, 1./freqr, 1./freqy);
      setup_f84gap(tblr, vr , alpha, gamma, 
		   freqa, freqc, freqg, freqt, freq_, 1./freqr, 1./freqy);
      for (i = 0; i < seq->endsite; i++)
        {
	  xx1 = &(q->x.s[i][0]);
	  xx2 = &(r->x.s[i][0]);
	  xx3 = &(mother->x.s[i][0]);
          
	  xx1t0 = (*xx1)[0];
	  xx1t1 = (*xx1)[1];
	  xx1t2 = (*xx1)[2];
	  xx1t3 = (*xx1)[3];
	  xx1t4 = (*xx1)[4];
          
	  xx2t0 = (*xx2)[0];
	  xx2t1 = (*xx2)[1];
	  xx2t2 = (*xx2)[2];
	  xx2t3 = (*xx2)[3];
	  xx2t4 = (*xx2)[4];
            
	  sa1 = tblq[0] * xx1t0 + tblq[1] * xx1t1 + 
	    tblq[2] * xx1t2 + tblq[3] * xx1t3 +
	    tblq[4] * xx1t4;
	  sa2 = tblr[0] * xx2t0 + tblr[1] * xx2t1 + 
	    tblr[1] * xx2t2 + tblr[3] * xx2t3 + 
	    tblr[4] * xx2t4;
	  (*xx3)[0] = sa1 * sa2;

            
	  sc1 = tblq[5] * xx1t0 + tblq[6] * xx1t1 + 
	    tblq[7] * xx1t2 + tblq[8] * xx1t3 +
	    tblq[9] * xx1t4;
	  sc2 = tblr[5] * xx2t0 + tblr[6] * xx2t1 + 
	    tblr[6] * xx2t2 + tblr[8] * xx2t3 + 
	    tblr[9] * xx2t4;
	  (*xx3)[1] = sc1 * sc2;

            
	  sg1 = tblq[10] * xx1t0 + tblq[11] * xx1t1 + 
	    tblq[12] * xx1t2 + tblq[13] * xx1t3 +
	    tblq[14] * xx1t4;
	  sg2 = tblr[10] * xx2t0 + tblr[11] * xx2t1 + 
	    tblr[11] * xx2t2 + tblr[13] * xx2t3 + 
	    tblr[14] * xx2t4;
	  (*xx3)[2] = sg1 * sg2;

            
	  st1 = tblq[15] * xx1t0 + tblq[16] * xx1t1 + 
	    tblq[17] * xx1t2 + tblq[18] * xx1t3 +
	    tblq[19] * xx1t4;
	  st2 = tblr[15] * xx2t0 + tblr[16] * xx2t1 + 
	    tblr[16] * xx2t2 + tblr[18] * xx2t3 + 
	    tblr[19] * xx2t4;
	  (*xx3)[3] = st1 * st2;

            
	  s_1 = tblq[20] * xx1t0 + tblq[21] * xx1t1 + 
	    tblq[22] * xx1t2 + tblq[23] * xx1t3 +
	    tblq[24] * xx1t4;
	  s_2 = tblr[20] * xx2t0 + tblr[21] * xx2t1 + 
	    tblr[21] * xx2t2 + tblr[23] * xx2t3 + 
	    tblr[24] * xx2t4;
	  (*xx3)[4] = s_1 * s_2;

      mother->scale[i] = sxx1[i] + sxx2[i];
#ifdef DEBUG
	    //	printf("[(%f,%f,%f,%f) %f]",(*xx3)[0],(*xx3)[1],(*xx3)[2],(*xx3)[3],mother->scale[i]);
#endif
        }
#ifdef DEBUG
	//   printf("id=%li\n",mother->id);
#endif
        count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < seq->endsite; i++)
            {
                xx3 = &(mother->x.s[i][0]);
                sxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                sxx3m = MAX (sxx3m, (*xx3)[2]);
                sxx3m = MAX (sxx3m, (*xx3)[3]);
                (*xx3)[0] /= sxx3m, (*xx3)[1] /= sxx3m, (*xx3)[2] /= sxx3m,
                    (*xx3)[3] /= sxx3m;
                mother->scale[i] += LOG (sxx3m);
            }
        }
    }
    else
    {
      // site rate variation
      for (i = 0; i < rcategs; i++)
	for (j = 0; j < categs; j++)
	  {
	    tblij = tbl[i][j];
	    tblij->ww2 = EXP (tblij->ratxi * lw2);
	    tblij->zz2 = EXP (tblij->ratxv * lw2);
	    tblij->ww2zz2 = tblij->ww2 * tblij->zz2;
	    tblij->vv2zz2 = (1.0 - tblij->ww2) * tblij->zz2;
	    tblr = &transprobr[i * 5 + j];
	    tblq = &transprobq[i * 5 + j];
	    setup_f84gap(tblq, tblij->rat * vq ,tblij->ratxi, tblij->ratxv, 
			 freqa, freqc, freqg, freqt, freq_, 1./freqr, 1./freqy);
	    setup_f84gap(tblr, tblij->rat * vr , tblij->ratxi, tblij->ratxv, 
			 freqa, freqc, freqg, freqt, freq_, 1./freqr, 1./freqy);
	  }
      for (i = 0; i < seq->endsite; i++)
	{
	  k = seq->category[seq->alias[i] - 1] - 1;
	  for (j = 0; j < rcategs; j++)
	    {

	      tbljk = tbl[j][k];
	      xx1 = &(q->x.s[i][j]);
	      xx2 = &(r->x.s[i][j]);
	      xx3 = &(mother->x.s[i][j]);
              
              
	      xx1t0 = (*xx1)[0];
	      xx1t1 = (*xx1)[1];
	      xx1t2 = (*xx1)[2];
	      xx1t3 = (*xx1)[3];
	      xx1t4 = (*xx1)[4];
              
	      xx2t0 = (*xx2)[0];
	      xx2t1 = (*xx2)[1];
	      xx2t2 = (*xx2)[2];
	      xx2t3 = (*xx2)[3];
	      xx2t4 = (*xx2)[4];

            
	  sa1 = tblq[0] * xx1t0 + tblq[1] * xx1t1 + 
	    tblq[2] * xx1t2 + tblq[3] * xx1t3 +
	    tblq[4] * xx1t4;
	  sa2 = tblr[0] * xx2t0 + tblr[1] * xx2t1 + 
	    tblr[1] * xx2t2 + tblr[3] * xx2t3 + 
	    tblr[4] * xx2t4;
	  (*xx3)[0] = sa1 * sa2;

            
	  sc1 = tblq[5] * xx1t0 + tblq[6] * xx1t1 + 
	    tblq[7] * xx1t2 + tblq[8] * xx1t3 +
	    tblq[9] * xx1t4;
	  sc2 = tblr[5] * xx2t0 + tblr[6] * xx2t1 + 
	    tblr[6] * xx2t2 + tblr[8] * xx2t3 + 
	    tblr[9] * xx2t4;
	  (*xx3)[1] = sc1 * sc2;

            
	  sg1 = tblq[10] * xx1t0 + tblq[11] * xx1t1 + 
	    tblq[12] * xx1t2 + tblq[13] * xx1t3 +
	    tblq[14] * xx1t4;
	  sg2 = tblr[10] * xx2t0 + tblr[11] * xx2t1 + 
	    tblr[11] * xx2t2 + tblr[13] * xx2t3 + 
	    tblr[14] * xx2t4;
	  (*xx3)[2] = sg1 * sg2;

            
	  st1 = tblq[15] * xx1t0 + tblq[16] * xx1t1 + 
	    tblq[17] * xx1t2 + tblq[18] * xx1t3 +
	    tblq[19] * xx1t4;
	  st2 = tblr[15] * xx2t0 + tblr[16] * xx2t1 + 
	    tblr[16] * xx2t2 + tblr[18] * xx2t3 + 
	    tblr[19] * xx2t4;
	  (*xx3)[3] = st1 * st2;

            
	  s_1 = tblq[20] * xx1t0 + tblq[21] * xx1t1 + 
	    tblq[22] * xx1t2 + tblq[23] * xx1t3 +
	    tblq[24] * xx1t4;
	  s_2 = tblr[20] * xx2t0 + tblr[21] * xx2t1 + 
	    tblr[21] * xx2t2 + tblr[23] * xx2t3 + 
	    tblr[24] * xx2t4;
	  (*xx3)[4] = s_1 * s_2;


	    }
	  mother->scale[i] = sxx1[i] + sxx2[i];
	}
      count++;
        if (count == SCALEINTERVAL)
        {
            count = 0;
            for (i = 0; i < seq->endsite; i++)
            {
                sxx3m = -MYREAL_MAX;
                for (j = 0; j < rcategs; j++)
                {
                    xx3 = &(mother->x.s[i][j]);
                    tempsxx3m = MAX ((*xx3)[0], (*xx3)[1]);
                    tempsxx3m = MAX (sxx3m, (*xx3)[2]);
                    tempsxx3m = MAX (sxx3m, (*xx3)[3]);
                    if (tempsxx3m > sxx3m)
                        sxx3m = tempsxx3m;
                }
                for (j = 0; j < rcategs; j++)
                {
                    xx3 = &(mother->x.s[i][j]);
                    (*xx3)[0] /= sxx3m, (*xx3)[1] /= sxx3m, (*xx3)[2] /= sxx3m,
                        (*xx3)[3] /= sxx3m;
                }
                mother->scale[i] += LOG (sxx3m);
            }
        }
    }
}    /* nuview */
#endif /*GAP*/
