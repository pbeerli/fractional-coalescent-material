// code to find the number cores on a computer
// [code from David Swofford, July 2007]
//
//
#include <stdio.h>

int countProcessors(void)
{

#if defined(CTL_HW) && defined(HW_NCPU)

  int mib[2], nprocessors;
  size_t len;

  mib[0] = CTL_HW;
  mib[1] = HW_NCPU;
  len = sizeof(nprocessors);
  if (sysctl(mib, 2, &nprocessors, &len, NULL, 0) != 0)
    nprocessors = 0;
  return nprocessors;
#    else
  return -1;
#    endif
}


int main(long *argc, char ** argv)
{
  int  nprocessors;

#if defined(_SC_NPROCESSORS_ONLN)
  nprocessors = sysconf(_SC_NPROCESSORS_ONLN);
#else
  nprocessors = countProcessors();
#endif

  printf("%i\n",nprocessors);
  return 0;
}
