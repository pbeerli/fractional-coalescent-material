#ifndef MIGRATE_DEFINITION /*definitions.h*/
#define MIGRATE_DEFINITION
/*-----------------------------------------------------------------
  Maximum likelihood estimation of migration rates 
  using coalescent trees

  Peter Beerli
  Genetics 357360
  University of Washington
  Seattle, WA 98195-7360, USA
  beerli@genetics.washington.edu
 
  With help of Joe Felsenstein (joe@agenetics.washington.edu),
  Mary Kuhner and Jon Yamato (mkkuhner@genetics.washington.edu)

  *----------------------------------------------------------------
  */
#ifndef VERSION
#define VERSION "alpha.2"
#endif
/*------------------------------------------------
  Adjusting the source to the system configuration
*/
 

#ifdef HAVE_CONFIG_H
#include "conf.h"
#endif


#ifdef __MWERKS__
#define MAC
#define HAVE_LGAMMA 0
#define HAVE_STRFTIME 1
#define HIGHBITS
#define LAMARC_MALLOC
//#include <ansi_prefix.mac.h> /*fixes time problems*/
//#include "unix.mac.h"
#endif /*MWERKS*/

#ifdef __WATCOMC__
#define QUICKC
#define WATCOM
#define DOS
#define HAVE_LGAMMA 0
#define HAVE_STRFTIME 1
#define HIGHBITS
#define LAMARC_MALLOC
#endif /*__WATCOM__*/

/* we have found no strftime() and replace the time
   output with blanks*/
#if HAVE_STRFTIME==0
#define NOTIME_FUNC
#endif
/* we have found no lgamma() function on the system
   so we use our own */
#if HAVE_LGAMMA==0
#undef HAVE_LGAMMA
#endif

#ifdef  GNUDOS
#define DJGPP
#define DOS
#endif /*GNUDOS*/

#ifdef __CMS_OPEN
#define CMS
#define EBCDIC true
#define SEEDFILE "seedfile data"
#define PARMFILE "parmfile data"
#define INFILE "infile data"
#define WEIGHTFILE "weightfile data"
#define CATFILE "catfile data"
#define OUTFILE "outfile data"
#define INTREE "intree data"
#define OUTTREE "outtree data"
#define MATHFILE "mathfile data"
#define SUMFILE "sumfile data"
#else
#define EBCDIC false
#define INFILE "infile"
#define WEIGHTFILE "weightfile"
#define CATFILE "catfile"
#define SEEDFILE "seedfile"
#define PARMFILE "parmfile"
#define OUTFILE "outfile"
#define TREEFILE "treefile"
#define UTREEFILE "usertree"
#define INTREE "intree"
#define OUTTREE "outtree"
#define MATHFILE "mathfile"
#define SUMFILE "sumfile"
#endif /*CMS_OPEN*/

#ifdef L_ctermid            /* try and detect for sysV or V7. */
#define SYSTEM_FIVE
#endif /*L_ctermid*/

#ifdef sequent
#define SYSTEM_FIVE
#endif /*sequent*/

#ifndef MAC
#ifndef SYSTEM_FIVE
# include<stdlib.h>
# if defined(_STDLIB_H_) || defined(_H_STDLIB) || defined(H_SCCSID) || defined(unix)
# define UNIX
# define MACHINE_TYPE "BSD Unix C"
# endif /* defined......*/
#endif /*SYSTEM_FIVE*/
#endif /*MAC*/

#ifdef __STDIO_LOADED
#define VMS
#define MACHINE_TYPE "VAX/VMS C"
#define printf vax_printf_is_broken
#define fprintf vax_fprintf_is_broken
void vax_printf_is_broken(const char *fmt,...);
void vax_fprintf_is_broken(FILE *fp,const char *fmt,...);
void vax_tweak_fmt(char *);
#endif /*__STDIO_LOADED*/

#ifdef _QC
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#include "graph.h"
#define DOS
#endif /*_QC*/

#ifdef _DOS_MODE
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS           /* DOS is  always defined if  on a dos machine */
#define MSC           /* MSC is defined for microsoft C              */
#endif /*DOS_MODE*/

#ifdef __MSDOS__      /* TURBO c compiler, ONLY (no other DOS C compilers) */
#define DOS
#define TURBOC
#include<stdlib.h>
#include<graphics.h>
#endif /*__MSDOS__*/

#ifdef DJGPP          /* DJ's gnu  C/C++ port */
#include<graphics.h>
#endif /*DJGPP*/

#ifndef MACHINE_TYPE
#define MACHINE_TYPE "ANSI C"
#endif

#ifdef DOS
#define MALLOCRETURN void 
#else
#define MALLOCRETURN void
#endif /*DOS*/

#ifdef VMS
#define signed /* signed doesn't exist in VMS */
#endif

/* default screen types */
#ifdef DOS
#define IBMCRT true
#define ANSICRT false
#else
#ifdef MAC
#define IBMCRT false 
#define ANSICRT false
#else
#define IBMCRT false 
#define ANSICRT true 
#endif /*MAC*/
#endif /*DOS*/

#ifdef DJGPP
#undef MALLOCRETURN
#define MALLOCRETURN void
#endif


/* includes: */
#ifdef UNIX
#include <strings.h>
#else
#include <string.h>
#endif /*UNIX*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#ifdef MAC
#include "mac_interface.h"
#endif /*MAC*/

#include <time.h>
#include <math.h>
#include <limits.h>
#include <signal.h>

#ifdef __MWERKS__
#undef DBL_EPSILON
#undef DBL_MAX
#include <float.h>
#endif /*__MWERKS__*/
#ifdef HAVE_MALLOCWRAP
#define LAMARC_MALLOC /*if we use this then we guard against allocation problems*/
/* redefinitions of calloc, malloc, realloc are in sighandler.h*/
#endif
#include "sighandler.h"
#define FClose(file) if (file) fclose(file) ; file=NULL

typedef unsigned char boolean;

#ifndef TRUE
#define TRUE    1
#endif

#ifndef FALSE
#define FALSE   0
#endif

#define SETBITS 32

#ifdef MAC
MALLOCRETURN    *mymalloc(long);
#else
MALLOCRETURN    *mymalloc();
#endif /*MAC*/


/*=================================================================
  SPECIFIC FOR MIGRATE

  *----------------------------------------------------------------
  + = works
  - = not yet implemented (needed?)

  */


/* number of k int k-allele model, the number  
   MUST really be bigger than the number of 
   observed alleles, much bigger does not harm	*/
#define MAXALLELES            93
#define LINESIZE            1024
#define DEFAULT_NMLENGTH      10 /* length of individual names*/
#define DEFAULT_ALLELENMLENGTH 6 /* length of allele names*/
#define DEFAULT_POPNMLENGTH  100 /* length of world names*/ 
#define NUMPOP 2
#define TWOHUNDRED 200
/* don't change below here ------------------------------------
   some other constants */

/* numerical borders and epsilons*/
#define SMALL_VALUE       10e-21
#define EPSILON         0.000001   /* a small number */
#define SMALLEPSILON       1e-15 /* a smaller number */
#define BIGEPSILON         0.001   /* a not so small number */
#define LONGCHAINEPSILON     100   /* an unreasonable big number, 
                                      so that the given number of 
                                      long chains is used*/
#ifndef MAXLONG
#define MAXLONG ((long)0x7fffffff)
#endif
#ifndef DBL_EPSILON
#include <float.h>
#ifndef DBL_EPSILON
#define DBL_MAX ((double)1.7976931348623157e308)
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif
#define PRECISION 0.0000001
/* some math constants */
#define LOGDBL_EPSILON -36.04365338911715608 /*N[Log[DBL_EPSILON] ,30]*/
#define LOG2 0.693147180559945309417232121458 /*N[Log[2] ,30]*/
#define LOG2PIHALF -0.918938533204672741780329736406 /*N[Log[1/Sqrt[2 Pi]] ,30]*/
#define TWOPI 6.28318530717958647692528676656 /*N[Log[2 Pi] ,30]*/
#define ROOTLENGTH         10000
#define FIRSTCHAIN            -1
#define FIRSTSTEP             -1
/* min/max functions */
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#define MIN3(a,b,c) (((a)<(b))&&((a)<(c)) ? (a) : (((b)<(c)) ? (b) : (c)))
#define MAX3(a,b,c) (((a)>(b))&&((a)>(c)) ? (a) : (((b)>(c)) ? (b) : (c)))


/* theta (4 Ne mu) related material */
#define SMALLEST_THETA	   1e-10
#define BIGGEST_THETA	    10e3
#define SMALLEST_MIGRATION 1e-10
#define BIGGEST_MIGRATION   10e9
#define SMALLEST_PROB     1e-100
#define SMALLEST_GAMMA      1e-3
#define BIGGEST_GAMMA        1e9
#define FRACTION_ALONG      0.66
#define TIMELIST_GUESS       100 
#define SAMPLETREE_GUESS       2 
#define START_ALPHA          1.0
/* random material */
#define AUTO                   0 /*+ use time() for seed*/
#define NOAUTO                 1 /*+ seed in parmfile   */
#define NOAUTOSELF             2 /*+ seed in seedfile   */

/* definitions for the first theta0 values */
#define OWN                    0 /*+ start values in parmfile*/
#define WATTERSON              1 /*- start values are Watterson estimate*/
#define EWENS                  2 /*- start values are Ewens estimate*/
#define FST                    3 /*+ start values are FST estimate*/
/* defines for the migration0 values */
/* define OWN                  0   + start values in parmfile */
/* define FST                  3   + start values are FST estimates*/
#define SLATKIN                1 /*- start values with slatkin's method*/
/* definitions for the type of FST we use*/
#define THETAVARIABLE         'T'   
#define MVARIABLE             'M'
/* definitions for the sankoff procedure */
/* the SANKOFF_DELTA should be smaller than any cost_ij*/
#define SANKOFF_DELTA 0.1

/* Newton-Raphson and gamma things */
#define NTRIALS             1000
#define LOCI_NORM          0.00001
#define GAMMA_INTERVALS      100

/* definitions of migration model */
#define MATRIX                 0 /*+ */    
#define ISLAND                 1 /*+ */    
#define STEPSTONE              2 /*- */
#define CONTINUUM              3 /*- */ 
#define NEIGHBOR               4 /*- */

/* print and plot material*/
#define PLOTALL                0 /*+ plot to outfile and mathfile*/
#define PLOTOUTFILE            1 /*+ plot only to outifle*/
/* more plotstuff in world.c */
#define PAGEFEED               fprintf(outfile,"\n\f\n")
#define PAGEFEEDWORLD               fprintf(world->outfile,"\n\f\n")
/* tree-print options */
#define NONE                   0 /*+ */
#define ALL                    1 /*+ */
#define BEST                   2 /*+ */
#define LASTCHAIN              3 /*+ */
/* likelihood ration stuff*/
#define LRATIO_STRINGS      1000 /*+*/

/* profile likelihood stuff*/
/* ALL and NONE are alread defined*/
#define TABLES                 2 /*+*/
#define SUMMARY                3 /*+*/
#define MAX_PROFILE_TRIALS   100 /*+*/
/* likelihood ratio stuff*/
#define MEAN 0
#define LOCUS 1
#define MAXPOP   10      /*n*(n-1)+n)*/
#define MAXPARAM 100
#define HEADER 1

/* stop while loops after some time*/
#define PANIC_MAX           1000

/* microsatellite stuff*/
#define MICRO_THRESHOLD       10
#define MAX_MICROSTEPNUM     200
#define XBROWN_SIZE            3
/* sequence stuff*/
#define MAXCATEGS              9
#define MANYCATEGS             2
#define ONECATEG               1
#endif /*definitions*/






