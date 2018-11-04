#include <migration.h>

void sortdata2(char ***result, data_fmt *data, world_fmt *world, long sublocus);
{
  for(pop=0; pop < data->numpop; pop++)
    {
      for(ind=0; ind < data->numind[pop][locus]; ind++)
	{
	  for(j=0; j < world->mutationmodel[sublocus].numsites; j++)
	    {
	      data->yy[pop][ind][locus][0][j],
		  }
  transpose_yy(result,data->yy,data, sublocus);
  //a = generate columnxrow matrix
  //qsort(a)
}


void sort_data(world_fmt *world, data_fmt *data, option_fmt *options, long locus)
{
  //foreach mutationmodel
  //   sortdata
  //      read from yy into simple nxm array elements are strings
  //      transpose data matrix
  //      sort
  //      transpose data matrix
  //   adjust and calculate weights
  long sublocus;
  for(sublocus=world->sublocistarts[locus]; sublocus < world->sublocistarts[locus]+data->subloci[locus];sublocus++)
    {
      sortdata2(data, world, sublocus);
      //adjust_weights(data, world, sublocus);
    }
}

void simplify(char **d, data_fmt *data, option_fmt *options, world_fmt *world, long sublocus, long locus)
{
  long z=0;
  long ii;
  long k;
  long pop;
  long zpop=0;
  long ind;
  long top;
  long numsites = world->mutationmodels[sublocus].numsites;
  enum siteclass_enum dataclass = world->mutationmodels[sublocus].dataclass;
  long datasize;
  
  if(dataclass==SITEWORD)
    {  
      datasize = options->allelenmlength;
    }
  else
    {
      datasize = 1;
    }

  world->sumtips = 0;

  for (pop = 0; pop < data->numpop; pop++)
    {
      zpop = 0;
      if(options->randomsubset>0 && options->randomsubset <  data->numind[pop][locus])
	{	
	  top = options->randomsubset;
	}
      else
	{
	  top = data->numind[pop][locus];
	}
      for (ii = 0; ii < top; ii++)
	{
	  ind = data->shuffled[pop][locus][ii];
	  zpop++;
	  for(k=0; k < numsites; k++)
	    {
	      memcpy(d[k][z], data->yy[pop][ind][locus][0][k],sizeof(char)*datasize);
	    }
	  z++;
        }
      data->numalleles[pop][locus] = zpop;
      world->sumtips += zpop;
    }
  qsort((void *) d, numsites, sizeof(char **), strcmp);
}

void sortdata2(data_fmt *data, world_fmt *world, long sublocus)
{
  char **d;
  long nsites = world->mutationmodels[sublocus].nsites;
  charvec2d(&d,nsites, world->sumtips);//this is the transpose
  simplify_data(&d, data, world);
}
