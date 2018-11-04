#include "migration.h"
typedef struct _custmmig {
  char * custm;    // *0d***00*
  char * custm2;   // ***0d**00  
  unsigned int  * popindex; // [ 0 1 2 0 0 1 1 2 2 ]
  unsigned int  * pindex; //   [ 0 1 2 3 4 5 6 7 8 ]
  unsigned int  * pclass; // 0=not estimated, 1=estimated, 2=symmetric, 4=+, 8=-, 16=estimated Nm,
                          //32=symmetric Nm, 64=+Nm, 128=-Nm, 256=divergence estimated, 
                          // 512=theta fixed,1024=M fixed, 2048=Nm fixed, 4096=Divergence fixed
  // examples: 0=not estimated,
  //           1=estimated 
  //           256+1=divergence+migration
  //           256+4=divergence + larger migration than other direction, this will often pair with
  //           256+8 or 256+1
  // 0=0
  // 1=1
  // 2=10
  // 4=100
  // 8=1000
  unsigned int  * pfunc;  //   0=uniform,1=exponential, 2=gamma
  unsigned int  * ptype;  // 1=theta, 2=migration, 4=divergence
  // this should be able to deliver an
  // TRUE in a comparison for example: is_mean, is_zero, is_smaller, is_bigger
  // can we have distribution of migration changes from speciation event towards today
  // or distribution from today into the past
  // for example instead of uniform migration events proposal these would suggest proposals 
  // according to the distribution, this is not not thought through yet.
} custm_migration_matrx_fmt;


boolean is_filler_custm(char ch)
{
  switch(ch)
    {
    case '[':
    case ']':
    case ',':
    case ' ':
    case '\t':
    case '\0':
      return TRUE;
    default:
      return FALSE;
    }
  return FALSE;
}

char is_content_custm(char ch)
{
  switch(ch)
    {
    case 'S':
    case 'M':
      return ch;
      break;
    default:
      return tolower(ch);
    }
  return ch;
}

char read_custom_matrix(FILE *file, option_fmt *options, char *value, long customnumpop, int *z)
{
  // using JASON
  // value is filled with 'custom-migration = { [.....]<,\n\r\t>{variable=keyword+keyword} {..} }'
  //find [
  int zz = 0;
  while (value[*z]!='\0' && value[*z]!='[')
    z++;
  if(value[*z]=='\0')
    {
      FGETS(value, LINESIZE, file);
      *z=0;
    }
  while (value[*z]!='\0' && value[*z]!=']')
    {
      char ch = value[*z];
      if(!is_filler_custm(ch))
	{
	  if (zz>=customnumpop)
	    {
	      options->custmalloc += TEN;
	      options->custm  = (char *) myrealloc (options->custm, sizeof (char) * options->custmalloc);
	    }
	  options->custm[zz] = is_content_custm(ch);
	}
      *z++;
    }
  return value[*z];
}

void set_custm_key(int * keys,int ii,char * tmp)
{
  /*
    uniform risk of event: 
     theta     1 0001
     migrate   2 0010
     split     4 0100
     join      8 1000
    possible combinations: theta+migrate, split+migrate, join+migrate
    extensions not coded, exponential risk today->past: migration;
        exponential risk past->today: split+migrate
  */
  char ch = temp[0];
  switch(ch)
    {
    case 't':
      keys[ii] += 1;
      break;
    case 'm':
      keys[ii] += 2;
      break;
    case 's':
      keys[ii] += 4;
      break;
    case 'j':
      keys[ii] += 8;
      break;
    default:
      error("Revisit the custom migration assignments, error in specification of parameter interaction");
    }
}

boolean read_custom_migration(FILE *file, option_fmt *options, char *value, long customnumpop)
{
  // using JASON
  // value is filled with 'custom-migration = { [.....]<,\n\r\t>{variable=keyword+keyword} {..} }'
  //find [
  int ii;
  int z  = 0;
  int zz = 0;
  char ch;
  while (value[z]!='\0' && value[z]!='[')
    z++;
  if(value[z]=='\0')
    {
      FGETS(value, LINESIZE, file);
      z=0;
    }
  ch = read_custom_matrix(file, options, value, customnumpop, &z);
  FGETS(value, LINESIZE, file);
  ii = 0;
  variables = (char *) calloc(options->custmalloc,sizeof(char));
  keys = (int *) calloc(options->custmalloc,sizeof(char));
  while(strcmp(value,'}'))
    {
      z=0;
      if (value[0]=='c' && value[1]=='u')    
	{
	  char *v = strdup(value);
	  char *vptr = v;
	  tmp = strtok(v," ");//the keyword: custom
	  tmp = strtok(v,"=");//the variable name
	  sprintf(variables[ii],"%s",tmp);
	  tmp = strtok(v," ");//the key name
	  set_custm_key(keys,ii,tmp);
	  iii=2;
	  while(v!=NULL)
	    {
	      tmp = strtok(v," ");//the second/third key name
	      set_custm_key(keys,ii,tmp);
	    }
	  free(vptr);
	}
      FGETS(value, LINESIZE, file);
    }
  /*
    variables=keys   example custom a=theta(0.01)+migration(estimate)   short custom a=t(0.01)+m(e)
                             custom b=migration(0|estimate)             short custom b=m(0|e)
                             custom c=divergence(estimate)+migration(+) short custom c=d(e)+m(+)
                             custom i=divergence(estimate)               
    migration-matrix         * a a a
                             b * c b
			     b i * 0
			     b 0 0 *
   */
}
