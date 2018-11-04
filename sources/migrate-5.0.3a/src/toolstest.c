#include "migration.h"
///
/// trim white space on both sides, keep pointer
/// example "   this is a test      "
/// find total number of characters:23, start of non-whitespace: 3 (start at zero)
/// end is at 13
void trim(char **line)
{
  char *tmp = *line;
  long start=0;
  const char whitespace[]=" \t\n\r";
  while (isspace(*tmp)) 
    {
      start++;
      tmp++;
    }
  long total = strlen(*line);
  long end = total;
  while (isspace(*(*line + end - 1))) 
    {
      end--;
    }
  memmove((*line),(*line)+start,end-start);
  (*line)[end-start]='\0';
}

int main()
{
  printf("testing trim()\n");
  char * line = (char *) calloc(100, sizeof(char));
  strcpy(line,"   this has trailing and leading blanks   ");
  printf("original   :|%s|\n",line);
  trim(&line);
  printf("manipulated:|%s|\n",line);
  strcpy(line,"   this has NO trailing and leading blanks");
  printf("original   :|%s|\n",line);
  trim(&line);
  printf("manipulated:|%s|\n",line);
  strcpy(line,"this has trailing and NO leading blanks   ");
  printf("original   :|%s|\n",line);
  trim(&line);
  printf("manipulated:|%s|\n",line);
  free(line);
}
