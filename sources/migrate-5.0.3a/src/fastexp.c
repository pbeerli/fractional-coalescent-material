#include "migration.h"
#include "tools.h"
#include <math.h>

//#define myEXPA (1048576 / 0.693147180559945309417232121458)
#define myEXPA 1512775.39519518569383584038231 
#define myEXPC 60801 
static union
    { 
        double d; 
      //LITTLE ENDIAN
        struct 
	{ 
	  int j, i; 
	} n; 
    } _eco;

double fast_exp(double y) 
{ 
  const double a = 1048576/M_LN2;
  const double b_c = 1072632447; /* 1072693248 - 60801 */
  _eco.n.i = (int)(a*y + b_c);
 
  return _eco.d;
}




int main(int argc, char **argv)
{
  long i,ii=0;
  double x;
  double y;
  long max=1000000000;
  for(i=0;i<max;i++)
    {
      if(ii > 10)
	ii=0;
      //x = fast_exp((double)-(ii++));
                  y = exp((double)-(ii++));
      //printf("%10li %10.10f %10.10f\n",-i, x, y);
    }
  return 0;
}

