/*
  mig-histogram code
 
  Peter Beerli
  beerli@fsu.edu
 
    Copyright 2002 Peter Beerli
 
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
 
  $Id: mig-histogram.c 2158 2013-04-29 01:56:20Z beerli $
 
 */
/*! \file mig-histogram.c
Calculation and printing of migration event histogram
*/

#include "tools.h"
#include "mig-histogram.h"
#include "sighandler.h"
#include "migrate_mpi.h"
#ifdef PRETTY
#include "pretty.h"
#endif

void setup_mighist (world_fmt * world, option_fmt * options);
void print_mighist (world_fmt * world);
void
setup_plotfield (plotfield_fmt * plotfield, char thisplotype,
                 long xsize, long ysize, char xlabel[], char ylabel[],
                 char yflabel[], char title[], boolean print);

long calc_migtable (MYREAL **migtable, histogram_fmt * histogram,
                    mighistloci_fmt * aa, world_fmt * world, long loci);
MYREAL average (MYREAL *vec, long size, long *weight, MYREAL *se, long *n);
MYREAL calc_quantile (MYREAL *vec, long size, long *weight, MYREAL quantile);

#define NBINS 30
void
print_histogram_ascii (FILE * out, histogram_fmt ** histogram,
                       plotfield_fmt ** plotfield, long loci, long nmigs,
                       long bins, long *sum, MYREAL ***migtable);

void
prepare_hist (histogram_fmt * hist, MYREAL *mtime, long count, long *weight);

void minmax (histogram_fmt * hist, float *tempmin, float *tempmax);


void
print_mighist_file (FILE * mighist, world_fmt * world)
{
    mighistloci_fmt *aa;
    //long          copies;
    long loc, j, i;
    for (loc = 0; loc < world->loci; loc++)
    {
        aa = &world->mighistloci[loc];
        for (j = 0; j < aa->mighistnum; j++)
        {
            //copies = aa->mighist[j].copies;
            for (i = 0; i < aa->mighist[j].migeventsize; i++)
            {
                FPRINTF (mighist, "%c %25.20f %i %i %li T%010li L%010li\n",
// abbreviation m or c 
                         aa->mighist[j].migevents[i].from == aa->mighist[j].migevents[i].to ? 'c' : 'm', /* 'd' may need to be included here*/
                         aa->mighist[j].migevents[i].age,    // time 
                         aa->mighist[j].migevents[i].from,    // from population 
                         aa->mighist[j].migevents[i].to,    // to population 
                         aa->mighist[j].migevents[i].sumlines,    // sumlines
                         j,                                 // tree number
                         loc                                // locus
                         );
            }
        }
    }
}

void
calc_mighistvalues (world_fmt * world, MYREAL ***migtable,
                    histogram_fmt ** histogram, long *total)
{
    mighistloci_fmt *aa;
    long loc, p1;
    long loci1 = world->loci == 1 ? 1 : world->loci + 1;
    for (loc = 0; loc < loci1; loc++)
    {
        for (p1 = 0; p1 < world->numpop2; p1++)
            migtable[loc][p1][2] = NOAVERAGE;
        aa = &world->mighistloci[loc];
        if (loc == world->loci)
        {
            if (world->loci != 1)
                calc_migtable (migtable[loc], histogram[loc], world->mighistloci,
                               world, world->loci);
            else
                break;
        }
        else
            total[loc] =
                calc_migtable (migtable[loc], histogram[loc], aa, world, 1);
    }
}

void
print_mighist_output (FILE * out, world_fmt * world, MYREAL ***migtable,
                      long *total)
{
    long loc, p1;
    long loci1 = world->loci == 1 ? 1 : world->loci + 1;
    FPRINTF (out, "\n\nSummary of %s events\n", 
             world->options->mighist_all ? "coalescence and migration" : "migration");
    FPRINTF (out, "===============%s================\n\n", 
             world->options->mighist_all ? "=========================" : "=========");
    for (loc = 0; loc < world->loci; loc++)
        total[world->loci] += total[loc];
    for (loc = 0; loc < loci1; loc++)
    {    /* Each locus + Summary */
        if (loc != world->loci)
            FPRINTF (out, "Locus %li\n", loc + 1);
        else
            FPRINTF (out, "Over all loci\n");
        FPRINTF (out,
                 "---------------------------------------------------------\n");
        FPRINTF (out,
                 "Population   Time                             Frequency\n");
        FPRINTF (out, "             -----------------------------\n");
        FPRINTF (out, "From    To   Average    Median     SE\n");
        FPRINTF (out,
                 "---------------------------------------------------------\n");
        for (p1 = (world->options->mighist_all ? 0 : world->numpop); p1 < world->numpop2; p1++)
        {
                if (migtable[loc][p1][2] == NOAVERAGE)
                {
                    FPRINTF (out, "%4li %4li    No %s event encountered\n",
			     (long) migtable[loc][p1][0] + 1,
			     (long) migtable[loc][p1][1] + 1,
			     (world->options->mighist_all && p1 < world->numpop) ?
		      "coalescence" : "migration");
                }
                else
                {
                    FPRINTF (out,
                             "%4li %4li    %3.5f    %3.5f    %3.5f    %3.5f\n",
                             (long) migtable[loc][p1][0] + 1,
                             (long) migtable[loc][p1][1] + 1,
                             migtable[loc][p1][2], migtable[loc][p1][3],
                             migtable[loc][p1][4],
                             migtable[loc][p1][5] / total[loc]);
                }
        }
        FPRINTF (out,
                 "---------------------------------------------------------\n");
        FPRINTF (out, "\n");
    }
}

/*
 * print_mighist() prints a table with the frequency of migrations events
 * from and to per timeinterval that is 1/100 of the full time  that goes
 * from zero to the Maximum time in the record.
 *
 * PopFrom PopTo  Average-Time Median-Time SE "Probability"
 *
 */
void
print_mighist (world_fmt * world)
{
    long loc, i;
    char plotype;
    char xlabel[255];
    char ylabel[255];
    char yflabel[255];
    char title[255];
    long xsize;
    long ysize;
    long to;
    long from;
    FILE *out = world->outfile;
    long loci1 = world->loci == 1 ? 1 : world->loci + 1;
    long numpop = world->numpop;
    long numpop2 = world->numpop2;
    long *total;
    MYREAL ***migtable;
    
    //loci x numpop2 x {mean, median, se}

    plotfield_fmt **plotfield;
    histogram_fmt **histogram;
    //only for overall loci
    if (world->options->mighist)
    {
        total = (long *) mycalloc (loci1 + 1, sizeof (long));
        plotfield = (plotfield_fmt **) mycalloc (loci1, sizeof (plotfield_fmt *));
        migtable = (MYREAL ***) mycalloc (loci1, sizeof (MYREAL **));
        histogram = (histogram_fmt **) mycalloc (loci1, sizeof (histogram_fmt *));
        for (loc = 0; loc < loci1; ++loc)
        {
            plotfield[loc] =
                (plotfield_fmt *) mycalloc (numpop2, sizeof (plotfield_fmt));
            migtable[loc] =
                (MYREAL **) mycalloc (numpop2, sizeof (MYREAL *));
            histogram[loc] =
                (histogram_fmt *) mycalloc (numpop2, sizeof (histogram_fmt));
            for (i = 0; i < world->numpop2; ++i)
            {
                histogram[loc][i].time = NULL;
                histogram[loc][i].weight = NULL;
            }
            for (i = 0; i < world->numpop2; ++i)
            {
                migtable[loc][i] = (MYREAL *) mycalloc (6, sizeof (MYREAL));
            }
        }

        //setup histogram
        plotype = 'a';
        xsize = NBINS;
        ysize = MIGHIST_YSIZE;
        strcpy (xlabel, "Time");
        strcpy (yflabel, "Frequency");
        strcpy (ylabel, "Count");
        for (loc = 0; loc < loci1; loc++)
        {
            for (to = 0; to < numpop; ++to)
            {
                if(world->options->mighist_all)
                {
                    sprintf (title, "Coalescences in population %li",
                             to + 1);
                    setup_plotfield (&plotfield[loc][mm2m(to,to,numpop)], plotype, xsize, ysize,
                                     xlabel, ylabel, yflabel, title, TRUE);
                }
                for (from = 0; from < numpop; ++from)
                {
                    if(from!=to)
                    {
                        sprintf (title, "Migrations from population %li to %li",from + 1, to + 1);
                        setup_plotfield (&plotfield[loc][mm2m(from,to,numpop)], plotype, xsize, ysize,
                                         xlabel, ylabel, yflabel, title, TRUE);
                    }
                }
            }
        }
        print_mighist_file (world->mighistfile, world);
        calc_mighistvalues (world, migtable, histogram, total);
        print_mighist_output (world->outfile, world, migtable, total);

        print_histogram_ascii (out, histogram, plotfield, loci1, numpop2,
                               NBINS, total, migtable);
#ifdef PRETTY
	pdf_print_mighist_table(world, migtable, total);
        pdf_mig_histogram(histogram, plotfield, loci1, numpop2, NBINS, total, migtable, TRUE, world);
#endif
        //fflush (out);
        myfree(total);
        for (loc = 0; loc < loci1; loc++)
        {
            for (i = (world->options->mighist_all ? 0 : world->numpop); i < numpop2; ++i)
            {
	      myfree(histogram[loc][i].time);
	      myfree(histogram[loc][i].weight);
	      myfree(migtable[loc][i]);
	      myfree(plotfield[loc][i].data[0]);
	      myfree(plotfield[loc][i].data);
	      myfree(plotfield[loc][i].y);
	      myfree(plotfield[loc][i].yfreq);
            }
            myfree(migtable[loc]);
            myfree(plotfield[loc]);
            myfree(histogram[loc]);
        }
        myfree(migtable);
        myfree(plotfield);
        myfree(histogram);
    }
}

void
setup_plotfield (plotfield_fmt * plotfield, char thisplotype,
                 long xsize, long ysize, char xlabel[], char ylabel[],
                 char yflabel[], char title[], boolean print)
{
    long i;
    plotfield->print = print;
    plotfield->type = thisplotype;
    plotfield->xsize = xsize;
    plotfield->ysize = ysize;
    sprintf (plotfield->xaxis,"%-254.254s", xlabel);
    sprintf (plotfield->yaxis,"%-254.254s", ylabel);
    sprintf (plotfield->yfaxis,"%-254.254s", yflabel);
    sprintf (plotfield->title,"%-254.254s", title);
    plotfield->yfreq = (float *) mycalloc (ysize, sizeof (float));
    plotfield->y = (long *) mycalloc (ysize, sizeof (long));
    plotfield->data = (char **) mymalloc (sizeof (char *) * xsize);
    plotfield->data[0] = (char *) mycalloc (xsize * (ysize + 1), sizeof (char));
    for (i = 1; i < xsize; i++)
    {
        plotfield->data[i] = plotfield->data[0] + i * (ysize + 1);
    }
}



long
calc_migtable (MYREAL **migtable, histogram_fmt * histogram,
               mighistloci_fmt * aa, world_fmt * world, long loci)
{
    long p1, p2, pa, i, j;
    long maxloci, locus;
    long copies;
    long n, total = 0, maxsize;
    MYREAL se;
    MYREAL ***migtime;
    long ***gencount;
    long **migcount;
    long **size;
    migtime = (MYREAL ***) mycalloc (world->numpop, sizeof (MYREAL **));
    gencount = (long ***) mycalloc (world->numpop, sizeof (long **));
    migcount = (long **) mycalloc (world->numpop, sizeof (long *));
    size = (long **) mycalloc (world->numpop, sizeof (long *));
    maxsize = 1;
    if (loci == 1)
        maxloci = 1;
    else
        maxloci = world->loci;
    for (locus = 0; locus < maxloci; locus++)
    {
        for (j = 0; j < aa[locus].mighistnum; j++)
        {
            if (maxsize < aa[locus].mighist[j].migeventsize)
                maxsize = aa[locus].mighist[j].migeventsize;
        }
    }
    for (p1 = 0; p1 < world->numpop; ++p1)
    {
        migtime[p1] = (MYREAL **) mycalloc (world->numpop, sizeof (MYREAL *));
        gencount[p1] = (long **) mycalloc (world->numpop, sizeof (long *));
        migcount[p1] = (long *) mycalloc (world->numpop, sizeof (long));
        size[p1] = (long *) mycalloc (world->numpop, sizeof (long));
        for (p2 = 0; p2 < world->numpop; ++p2)
        {
            migtime[p1][p2] = (MYREAL *) mycalloc (maxsize, sizeof (MYREAL));
            gencount[p1][p2] = (long *) mycalloc (maxsize, sizeof (long));
            size[p1][p2] = maxsize;
        }
    }
    for (locus = 0; locus < maxloci; locus++)
    {
        for (j = 0; j < aa[locus].mighistnum; j++)
        {
            copies = aa[locus].mighist[j].copies;
            for (i = 0; i < aa[locus].mighist[j].migeventsize; i++)
            {
                p1 = (long) aa[locus].mighist[j].migevents[i].from;
                p2 = (long) aa[locus].mighist[j].migevents[i].to;
                if (migcount[p1][p2] >= size[p1][p2])
                {
                    size[p1][p2] += 10;
                    gencount[p1][p2] =
                        (long *) myrealloc (gencount[p1][p2],
                                          sizeof (long) * size[p1][p2]);
                    memset (gencount[p1][p2] + migcount[p1][p2], 0,
                            sizeof (long) * 10);
                    migtime[p1][p2] =
                        (MYREAL *) myrealloc (migtime[p1][p2],
                                            sizeof (MYREAL) * size[p1][p2]);
                    memset (migtime[p1][p2] + migcount[p1][p2], 0,
                            sizeof (MYREAL) * 10);
                }
                gencount[p1][p2][migcount[p1][p2]] += copies;
                migtime[p1][p2][migcount[p1][p2]] +=
                    aa[locus].mighist[j].migevents[i].age;
                migcount[p1][p2] += 1;
            }
        }
    }
    for (p1 = 0; p1 < world->numpop; p1++)
    {
        if(world->options->mighist_all)
        {            
            migtable[p1][0] = p1;
            migtable[p1][1] = p1;
            migtable[p1][2] = average (migtime[p1][p1],
                                       migcount[p1][p1],
                                       gencount[p1][p1], &se, &n);
            migtable[p1][4] = se;
            migtable[p1][5] = n;
            migtable[p1][3] = calc_quantile (migtime[p1][p1],
                                             migcount[p1][p1],
                                             gencount[p1][p1], 0.5);
            prepare_hist (&histogram[p1], migtime[p1][p1], migcount[p1][p1],
                          gencount[p1][p1]);
        }
        for (p2 = 0; p2 < world->numpop; p2++)
        {
            pa = mm2m (p1, p2, world->numpop);
            migtable[pa][0] = p1;
            migtable[pa][1] = p2;
            migtable[pa][2] = average (migtime[p1][p2],
                                       migcount[p1][p2],
                                       gencount[p1][p2], &se, &n);
            migtable[pa][4] = se;
            migtable[pa][5] = n;
            migtable[pa][3] = calc_quantile (migtime[p1][p2],
                                        migcount[p1][p2],
                                        gencount[p1][p2], 0.5);
            prepare_hist (&histogram[pa], migtime[p1][p2], migcount[p1][p2],
                          gencount[p1][p2]);
        }
    }
    for (p1 = (world->options->mighist_all ? 0 : world->numpop); p1 < world->numpop2; p1++)
    {
        total += (long) migtable[p1][5];
    }
    for (p1 = 0; p1 < world->numpop; ++p1)
    {
        myfree(migcount[p1]);
        myfree(size[p1]);
        for (p2 = 0; p2 < world->numpop; ++p2)
        {
            myfree(migtime[p1][p2]);
            myfree(gencount[p1][p2]);
        }
        myfree(migtime[p1]);
        myfree(gencount[p1]);
    }
    myfree(migtime);
    myfree(gencount);
    return total;
}

void
prepare_hist (histogram_fmt * hist, MYREAL *mtime, long count, long *weight)
{
    hist->count = count;
    if (hist->time == NULL)
        hist->time = (MYREAL *) mycalloc (count + 1, sizeof (MYREAL));
    else
        hist->time = myrealloc (hist->time, sizeof (MYREAL) * (count + 1));
    if (hist->weight == NULL)
        hist->weight = (long *) mycalloc (count + 1, sizeof (long));
    else
        hist->weight = myrealloc (hist->weight, sizeof (long) * (count + 1));
    memcpy (hist->time, mtime, sizeof (MYREAL) * count);
    memcpy (hist->weight, weight, sizeof (long) * count);
}

MYREAL
average (MYREAL *vec, long size, long *weight, MYREAL *se, long *n)
{
    long i;
    MYREAL mean, sum = 0., sum2 = 0.;
    long sumweight = 0;
    for (i = 0; i < size; ++i)
        sumweight += weight[i];
    for (i = 0; i < size; ++i)
    {
        sum += vec[i] * weight[i];
        sum2 += (vec[i] * weight[i]) * (vec[i] * weight[i]);
    }
    if (sumweight != 0)
    {
        mean = sum / sumweight;
        if (sumweight > 1)
            *se = sqrt (fabs (sum - sum2)) / (sumweight - 1.);
        else
            *se = MYREAL_MAX;
        *n = sumweight;
        return mean;
    }
    else
    {
        *n = 0;
        *se = MYREAL_MAX;
    }
    return NOAVERAGE;
}

MYREAL
calc_quantile (MYREAL *vec, long size, long *weight, MYREAL quantile)
{
    long i, j, z = 0;
    MYREAL *tmp1;
    MYREAL val;
    long sumweight = 0;
    for (i = 0; i < size; ++i)
        sumweight += weight[i];

    tmp1 = (MYREAL *) mycalloc (sumweight + 1, sizeof (MYREAL));

    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < weight[i]; ++j)
            tmp1[z++] = vec[i];
    }
    qsort ((void *) tmp1, sumweight, sizeof (MYREAL), numcmp);
    val = tmp1[(long) (sumweight * quantile)];
    myfree(tmp1);
    return val;
}

void
increase_mighist (mighistloci_fmt * mighistlocus)
{
    long i;
    if (mighistlocus->allocsize <= mighistlocus->mighistnum + 1)
    {
        mighistlocus->allocsize += DEFAULTALLOCSIZE;
        mighistlocus->mighist = (mighist_fmt *)
                                 myrealloc (mighistlocus->mighist,
                                         sizeof (mighist_fmt) * mighistlocus->allocsize);
        for (i = mighistlocus->allocsize - DEFAULTALLOCSIZE;
                i < mighistlocus->allocsize; i++)
        {
            mighistlocus->mighist[i].migeventsize = 0;
	    mighistlocus->mighist[i].copies = 0;
	    mighistlocus->mighist[i].weight = 0;
            mighistlocus->mighist[i].allocsize = DEFAULTALLOCSIZE;
            mighistlocus->mighist[i].migevents =
                (migevent_fmt *) mycalloc (DEFAULTALLOCSIZE, sizeof (migevent_fmt));
        }
    }
}

void
setup_mighist (world_fmt * world, option_fmt * options)
{
    long locus, i;
    long allocsize = DEFAULTALLOCSIZE;
    if (world->options->mighist)
    {
        world->mighistloci = (mighistloci_fmt *)
                             mycalloc (world->loci+1, sizeof (mighistloci_fmt));
        world->mighistlocinum = 0;
        for (locus = 0; locus < world->loci+1; locus++)
        {
            world->mighistloci[locus].allocsize = allocsize;
            world->mighistloci[locus].mighistnum = 0;
            world->mighistloci[locus].mighist =
                (mighist_fmt *) mycalloc (world->mighistloci[locus].allocsize,
                                        sizeof (mighist_fmt));
            for (i = 0; i < allocsize; i++)
            {
                world->mighistloci[locus].mighist[i].migeventsize = 0;
                world->mighistloci[locus].mighist[i].allocsize = allocsize;
                world->mighistloci[locus].mighist[i].migevents =
                    (migevent_fmt *) mycalloc (allocsize, sizeof (migevent_fmt));
            }
        }
    }
    else
      world->mighistloci = NULL;
}

void
destroy_mighist (world_fmt * world)
{
    long locus, i;
    if (world->options->mighist)
    {
      for (locus = 0; locus < world->loci; locus++)
        {
	  
	  for (i = 0; i < world->mighistloci[locus].allocsize; i++)
            {
	      myfree(world->mighistloci[locus].mighist[i].migevents);
	      //	      myfree(world->mighistloci[locus].mighist[i].time);
            }
	  myfree(world->mighistloci[locus].mighist);
        }
      myfree(world->mighistloci);
    }
}

void
minmax (histogram_fmt * hist, float *tempmin, float *tempmax)
{
    long i;
    MYREAL tmp1, tmp2;
    MYREAL tmpmin = MYREAL_MAX;
    MYREAL tmpmax = -MYREAL_MAX;

    for (i = 0; i < hist->count; i++)
    {
        if ((tmp1 = hist->time[i]) < tmpmin)
            tmpmin = tmp1;
        if ((tmp2 = hist->time[i]) > tmpmax)
            tmpmax = tmp2;
    }
    *tempmax = (float) tmpmax;
    *tempmin = (float) tmpmin;
}


void
print_histogram_ascii (FILE * out, histogram_fmt ** histogram,
                       plotfield_fmt ** plotfield, long loci, long nmigs,
                       long bins, long *sum, MYREAL ***migtable)
{
    long loc, i, j, z, zz;
    float biggest = 0.0F;
    float *binning;
    float *binvec;
    float tempmin = FLT_MAX;
    float tempmax = -FLT_MAX;
    float begin = FLT_MAX;
    float end = -FLT_MAX;
    float delta;
    //MYREAL        sum = 0;
    MYREAL mtime;
    long weight;
    binning = (float *) mycalloc (bins+1, sizeof (float));
    binvec = (float *) mycalloc (bins+1, sizeof (float));

    for (loc = 0; loc < loci; loc++)
    {
        for (i = 0; i < nmigs; i++)
        {
            if (migtable[loc][i][2] == NOAVERAGE)
            {
                plotfield[loc][i].print = FALSE;
                continue;
                //no event for this migration from i to j
            }
            minmax (&histogram[loc][i], &tempmin, &tempmax);
            if (tempmin < begin)
                begin = tempmin;
            if (tempmax > end)
                end = tempmax;
        }
    }
    delta = (end - begin) / bins;
    binning[0] = begin + 0.5 * delta;
    for (i = 1; i < bins; i++)
        binning[i] = delta + binning[i - 1];
    for (loc = 0; loc < loci; loc++)
    {
        for (i = 0; i < nmigs; i++)
        {
            if ((migtable[loc][i][2] == NOAVERAGE) || (plotfield[loc][i].y == NULL))
                continue;
            //no event for this migration i->j
            memset (binvec, 0, sizeof (float) * (bins+1));
            for (j = 0; j < histogram[loc][i].count; j++)
            {
                mtime = histogram[loc][i].time[j];
                weight = histogram[loc][i].weight[j];
                z = 0;
                while (mtime > binning[z] && z < bins)
                    z++;
                binvec[z] += weight;
            }
            biggest = 0.;
            for (j = 0; j < bins; j++)
            {
                plotfield[loc][i].y[j] = (long) binvec[j];
                plotfield[loc][i].yfreq[j] = binvec[j] = binvec[j] / sum[loc];
                if (biggest < binvec[j])
                    biggest = binvec[j];
            }
            for (j = 0; j < bins; j++)
            {
                for (zz = 0;
                        zz <
                        (long) (binvec[j] * plotfield[loc][i].ysize / biggest);
                        zz++)
                    plotfield[loc][i].data[j][zz] = '+';
                plotfield[loc][i].data[j][zz] = '\0';
            }
        }
    }
    for (loc = 0; loc < loci; loc++)
    {
        if (loc == (loci - 1))
        {
            if (loci > 1)
                FPRINTF (out,
                         "\nOver all loci\n------------------------------------------------------------------\n");
            else
                FPRINTF (out,
                         "\nLocus %li\n------------------------------------------------------------------\n",
                         loc + 1);
        }
        else
            FPRINTF (out,
                     "\nLocus %li\n------------------------------------------------------------------\n",
                     loc + 1);

        for (i = 0; i < nmigs; i++)
        {
            if (plotfield[loc][i].print)
            {
                FPRINTF (out, "%s\n\n%10.10s %10.10s %10.10s\n",
                         plotfield[loc][i].title, plotfield[loc][i].xaxis,
                         plotfield[loc][i].yaxis, plotfield[loc][i].yfaxis);
                for (j = 0; j < bins; j++)
                {
		  /*                    FPRINTF (stdout,
                             "loc=%li i=%li j=%li %10.6f %10li %10.6f %s\n",
                             loc, i, j, binning[j], plotfield[loc][i].y[j],
                             plotfield[loc][i].yfreq[j],
                             plotfield[loc][i].data[j]);
		  */
                    FPRINTF (out,
                             "%10.6f %10li %10.6f %s\n", binning[j],
                             plotfield[loc][i].y[j], plotfield[loc][i].yfreq[j],
                             plotfield[loc][i].data[j]);
                }
                FPRINTF (out, " \n");
            }
        }
    }
}
