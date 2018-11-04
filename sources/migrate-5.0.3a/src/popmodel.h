#ifndef POPMODEL_INCLUDE
#define POPMODEL_INCLUDE
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
#include "migration.h"

///
/// enumeration of all poossible popmodel types
enum popmodel_enum
  {
    WRIGHTFISHER, CANNING, MORAN
  }

extern char[]          which_popmodel(option_fmt *options, char temp[]);



#endif /*POPMODEL_INCLUDE*/
