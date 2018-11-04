/* ! \file popmodel.c */
/*------------------------------------------------------
  Estimation 
  of migration rate  and effectice population size
  using a Metropolis-Hastings Monte Carlo algorithm
  -------------------------------------------------------
  P O P M O D E L    R O U T I N E S

  related to the use of different population models
  Peter Beerli 2006, Tallahassee
  beerli@fsu.edu

  Copyright 2006 Peter Beerli, Tallahassee FL

  This software is distributed free of charge for non-commercial use
  and is copyrighted. Of course, we do not guarantee that the software
  works and are not responsible for any damage you may cause or have.

  $Id$

  -------------------------------------------------------*/
#ifdef POPMODEL
#include "migration.h"
#include "popmodel.h"
#include "sighandler.h"
#include "tools.h"
#include "sequence.h"

#include "options.h"
#include "migrate_mpi.h"

#ifdef DMALLOC_FUNC_CHECK
#include <dmalloc.h>
#endif
/* prototypes ------------------------------------------- */
char[]          which_popmodel(option_fmt *options, char temp[]);


///
/// prints out the population model in use with its parameter
/// if any. Currently only three models are allowed. The function returns
/// the string with the name of the model in the parameter temp is also returned.
/// the string is_set is filled with menu information displayed in the popmodel menu 
char[]          which_popmodel(option_fmt *options, char temp[], char is_set[])
{
  switch(options->popmodel)
    {
    case MORAN:
      sprintf(temp,"Moran model");
      sprintf(is_set,"Set");
      break;
    case CANNING:
      sprintf(temp,"Canning model with sigma=%f",options->popmodel_sigma);
      sprintf(is_set,"%f",options->popmodel_sigma);
      break;
    case WRIGHTFISHER:
    default:
      sprintf(is_set,"Set");
      sprintf(temp,"Wright-Fisher model");
      break;
    }
  return temp;
}

///
/// set the population model from a list of possible models, currently
/// 3 models: Wright-Fisher (default), Canning with an additional parameter
/// the variance of the offspring number (Wright-Fisher is a special case of
/// the Canning model with variance=1), Moran model.
void set_popmodel(option_fmt *options)
{
  char is_set[LINESIZE];
  char temp[LINESIZE];
  which_popmodel(options,temp, is_set);
  do
    {
      char input[LINESIZE];
      fprintf(stdout,"Select a population model:\n");
      fprintf(stdout,"W     Wright-Fisher model          %4.4s\n",is_set);
      fprintf(stdout,"C     Canning's model              %4.4s\n",is_set);
      fprintf(stdout,"M     Moran's model                %4.4s\n",is_set);
      fprintf(stdout,"Q     return to previous menu\n\n");
      fprintf(stdout,"Type either W, C, M, or Q\n\n");
      FGETS(input, LINESIZE, stdin);
      if (input[0] == '\0')
	options->popmodel = WRIGHTFISHER;
      else
	{
	  switch(uppercase(input[0]))
	    {
	    case 'M':
	      sprintf(temp,"Moran model");
	      break;
	    case 'C':
	      fprintf(stdout,"What is the standard deviation (std) of the number of offspring?\n");
	      fprintf(stdout,"[Wright-Fisher has std=1, Models close to Moran's close to 1/N_e]\n> ");
	      
	      FGETS(input,LINESIZE,stdin);
	      if (input[0] == '\0')
		options->popmodel_sigma = 1.0;
	      else
		options->popmodel_sigma = atof(input);
	      
	      sprintf(temp,"Canning model with sigma=%f",options->popmodel_sigma);
	      break;
	    case 'Q':
	      input[0]='\0';
	      return;
	    case 'W':
	    default:
	      sprintf(temp,"Wright-Fisher model");
	      break;
	    }
	  input[0]='\0';
	}
    } while(1);
}
#endif /*POPMODEL*/
