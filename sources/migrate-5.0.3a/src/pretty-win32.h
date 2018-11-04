#ifndef PRETTY_WIN32
#define PRETTY_WIN32
#ifdef WINDOWS
#ifndef MPI
extern "C" {
  void set_haru_handler(void);
}
#else
extern void set_haru_handler(void);
#endif
#endif
#endif
