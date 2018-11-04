///////////////
// on windows systems this file is need when compiling migrate
// using the libharu PDF-package.
//
// Michal Palzcewski and Peter Beerli 2005
#ifdef WINDOWS
#include <new.h>
#ifndef MPI
extern "C" {
  void set_haru_handler(void);
}
extern int throw_new_handler(size_t t);
void set_haru_handler(void) 
{
  _set_new_handler(throw_new_handler);
}
#endif
#endif
int dummy_function(int do_nothing);

int dummy_function(int do_nothing)
{
  int x;
  int y = 1.0;
  //Some computer fail to compile MIGRATE with empty files
  x = do_nothing + y;
  return x;
}
