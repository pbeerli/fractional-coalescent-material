// generic table generator package
// containing functions:
// pdf_table()
// find_col_width()
// define_col_start()
// pdf_print_table_header()
 
///
/// finds the width of all columns given the elements of the table
void  find_col_width(int cols, int cols, char ***elements, char **header, float *col_widths)
{
  float w;
  float keepw;
  int row;
  int col;

  for(col=0;col < cols, col++)
    {
      keepw = (float) pdf_contents_get_text_width(canvas, header[col], NULL, NULL);
      for(row=0;row < rows; row++)
	{
	  w = (float) pdf_contents_get_text_width(canvas, elements[row][col], NULL, NULL);
	  if(w > keepw)
	    keepw = w;
	}
      col_width[col] = keepw;
    }
}

///
/// align column
float align_column(char position, float cw, float left_margin)
{
  float loc = 0.;
  switch(position)
    {
    case 'R':
      loc = left_margin + cw;
      break;
    case 'C':
      loc = left_margin + cw/2.;
      break;
    case 'L':
    default:
      loc = left_margin;
      break;
    }
  return loc;
}

///
/// sets X-coordinates to start the columns 
void  define_col_start(float cols, float * col_widths, int col_overflow, char *position, float left_margin, float page_width, float separator, float *col_starts)
{
  int col;

  col_starts[0] = align_column(position[0], col_widths[0], left_margin);
  for(col=1; col < cols; col++)
    {
      col_starts[col] =  align_column(position[col], col_width[col], separator + col_starts[col-1]);
      if((col_starts[col] + col_width[col]) > page_width) //correct for 'L', but questionable for 'R'
	{
	  if(col_overflow==0)
	    col_starts[col] = col_starts[0];
	  else
	    {
	      if(col_overflow < col)
		col_starts[col] = col_starts[col_overflow];
	      else
		colstarts[col] = col_starts[0];
	    }
	}
    }
}  

///
/// prints the table header at column positions
void  pdf_print_table_header(float *page_height, char * position, int cols, float *col_starts,char **header)
{
  int col;

  for(col=0; col < cols; col++)
    {
      pdf_printf(col_starts[col], *page_height, position[col], "%s", header[col]);
    }
}

///
/// generic table generator
/// \param float left_margin left edge of table
/// \param float * page_height page height 
/// \param int cols number of columns in table if there are too many columns then new line
/// \param int rows number of rows in table, if a new page is needed then header is repeated
/// \param char *** elements all table elements, currently no formatting of these
/// \param char ** header   header row
/// \param int col_overflow when there are too many columns this is the column to restart 
void pdf_table(float left_margin, float *page_height, int cols, int rows, char ***elements, char **header, int col_overflow)
{
  boolean new_page=FALSE;

  int row;
  int col;
  int *col_widths;
  
  float *col_starts;
  float oldcol = HUGE;

  col_widths = mycalloc(cols,sizeof(int));
  col_starts = mycalloc(cols,sizeof(float));

  find_col_width(cols, rows, elements, header, col_widths);
  define_col_start(col_widths,page_width);
  
  pdf_print_table_header(col_starts,header);
  pdf_draw_line(55,*page_height,right_margin, *page_height);
  for(row=0; row< rows; row++)
    {
      new_page = pdf_advance(page_height);
      if(new_page)
	{
	  pdf_print_table_header(col_starts,header);
	  pdf_draw_line(55,*page_height,right_margin, *page_height);
	  new_page = pdf_advance(page_height);
	}
      for(col=0;col < cols; col++)
	{
	  if(col_starts[col] < oldcol)
	    {
	      new_page = pdf_advance(page_height);
	      if(new_page)
		{
		  pdf_print_table_header(col_starts,header);
		  pdf_draw_line(55,*page_height,right_margin, *page_height);
		  new_page = pdf_advance(page_height);
		}
	    }
	  pdf_printf(col_starts[col],*page_height, position[col], "%s",elements[row][col]);
	  oldcol = col_starts[col];
	}
    }
  pdf_advance(page_height);
  myfree(col_widths);
  myfree(col_starts);
}
