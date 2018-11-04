#ifndef MIGRATION_HEADER
#define MIGRATION_HEADER

#include "definitions.h"
/*-----------------------------------------------------------------
  Maximum likelihood estimation of migration rates 
  using coalescent trees

  Peter Beerli
  Genetics 357360
  University of Washington
  Seattle, WA 98195-7360, USA
  beerli@genetics.washington.edu
 
  With help of Joe Felsenstein (joe@agenetics.washington.edu),
  Mary Kuhner and Jon Jamato (mkkuhner@genetics.washington.edu)

  *----------------------------------------------------------------
 
  Time, 'tyme', in this tree is measured from the tips to the root.
  I.e. the tips are at tyme '0', and the root node has the largest
  value for 'tyme'. 
  
  *----------------------------------------------------------------
  */

/* typedefs -----------------------------------------------------------*/
typedef long 		longer[6]; /* used in random.c */


/* defines the data structure read from infile*/
typedef struct _data {
	FILE *infile;
	char datatype;
	char *****data;
    long *maxalleles;
	char **popnames;
	long *numind;
	char ***indnames;
	long numpop;
	long loci;
	char dlm;
} data_fmt;

typedef struct _option {
	FILE *parmfile;
	FILE *seedfile;
    int menu;
    char dlm;
    long popnmlength;
    long nmlength;
    long allelenmlength;
    long autoseed;
    boolean interleaved;
    boolean printdata;
    boolean progress;
    boolean treeprint;
    double ttratio;
    double watt;
    boolean usertree;
    long ssteps;
    long sincrement;
    long schains;
    long lsteps;
    long lincrement;
    long lchains;
    double cxy[4];
    long thetastepping;
    long output_limit;
    double interp_cutoff;
    char datatype;
    long sites;
    int mc_method;
    boolean curve;
    boolean curve_tally;
    long curve_val;
    long curve_method;
    char curve_func;
	boolean autocorr;
	double alamda; /* for autocorrelation of DNA*/
} option_fmt;
	

/* used in the tree structure*/
typedef union  xarray_fmt {
	double *a;
	double *s;
	} xarray_fmt;

typedef struct _node {
  struct _node *next, *back;
  boolean tip;
  double probxi;
  double probxv;
  char type;
  long number;
  long pop;
  long actualpop;
  long id;
  xarray_fmt x;
  double lxmax;
  char *nayme;
  boolean top;
  boolean dirty;
  double v, tyme, length, xcoord;
  short ycoord, ymin, ymax;
  int horizonstate;
  long * mut_accum;
} node;

typedef struct vtlist {
    node *eventnode;   /* node with age=tyme */
    double age;        /* tyme from top nodelet */
    double interval;   /* interval t[i+1].age - t[i].age */
    long lineages[2];
    long from;
    long to;
    long pop;
    long slice;
} vtlist;

typedef struct _vtlist_array {
    long copies;
    long allocT;
    long T;
    long oldT;
    vtlist *tl;
} vtlist_array_fmt;

typedef struct tree_fmt {
  node   **nodep;
  node *root;
  long pop;
  long tips;
} tree_fmt ;        

typedef struct _world {
    tree_fmt * curtree;
    long *maxalleles;
    node *root;
    double *likelihood;
    double param_like;
    long locus;
    long numpop;
    long sumtips;
    long lineages[2];
    double *param0;
	vtlist_array_fmt *treetimes;
	FILE *outfile;
	FILE *curvefile;
	FILE *curvefile1;
	FILE *curvefile2;
	FILE *treefile;
	FILE *errfile;
} world_fmt;


typedef struct _nr_fmt {
    long partsize; /*number of part-variables, fixed per model*/
    double *parts; /* parts of the first and second derivatives*/
    double *d;     /* first derivates*/
    double **dd;   /* second derivatives*/
    double *param; /* changed values of param*/
    double *oparam;/* saved old parameters*/
    double llike;  /* parameter LOGlikelihood*/
    double ollike; /* old "   "*/
    double *datalike;/*P(D|G)*/
    double *apg0;  /* part-loglikelihoods of param0*/
    double *apg;   /* part-loglikelihoods of param*/
    double apg_max;/* maxvalue of apg*/
    long numpop;   /* number of worlds*/
    long numpop2;  /* 2*numpop*/
    double PGC;    /* uncorrected llike*/
    double oPGC;    /* uncorrected ollike*/
} nr_fmt;


typedef struct _migr_table_fmt {
    long from;
    long to;
    double time;
} migr_table_fmt;   

typedef struct proposal_fmt {
    node *origin;
    node *target; 
    node *realtarget; 
    node *tsister;
    node *realtsister;
    node *osister;
    node *realosister;
    node *ocousin;
    node *realocousin;
    node *oback;
    node *realoback;
    node **line_f;
    node **line_t;
    node *connect;
    double likelihood;
    double time;
    double v;
    double vs;
    double *xf;
    double *xt;
    node **aboveorigin;
    node **bordernodes;
    migr_table_fmt *migr_table;
    migr_table_fmt *migr_table2;
    long migr_table_counter;
    long migr_table_counter2;
    long old_migr_table_counter;
    long old_migr_table_counter2;
    long timeslice;
} proposal_fmt;


/* global variables (this should be not used at all)*/
char application[LINESIZE];  
#endif








