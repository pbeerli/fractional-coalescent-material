void data_menu(data_fmt *data, option_fmt *options)
{
  boolean done= FALSE;
  char *input = mycalloc(LINESIZE,sizeof(char));
  while(!done)
    {
      get_data_menu_title(data, options);
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
