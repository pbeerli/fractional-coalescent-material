#include "migration.h"

void data_menu(data_fmt *data, option_fmt *options)
{
  boolean done= FALSE;
  char *input = mycalloc(LINESIZE,sizeof(char));
  while(!done)
    {
      get_menu_title("DATA MENU",data, options);
      get_datafile(data, options);
      get_auxdatafiles(data, options);
      get_loci(data, options);
      get_locuslist(data, options);
      get_mutationmodel(data, options);
      
      done = read_input(input);
      
      set_datafile(input,data, options);
      set_auxdatafiles(input, data, options);
      set_loci(input, data, options);
      set_locuslist(input, data, options);
      set_mutationmodel(input, data, options);
    }
  myfree(input);
}

void underline(long n)
{
  long i;
  for(i=0;i<n;i++)
    fprintf(stdout,"-");
  fprintf(stdout,"\n");
}

void get_data_menu_title(char title[], option_fmt *options, data_fmt *data)
{
  fprintf(stdout,"%s\n",title);
  underline("-",strlen(title));
}

void get_datafile(data_fmt * data, option_fmt *options)
{
  char * datafilestring = options->infilename;
  fprintf(stdout,"%li Datafile:%60.60s", 1, datafilestring);
}
