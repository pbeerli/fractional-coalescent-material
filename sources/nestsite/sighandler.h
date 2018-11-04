


#ifndef SIGNAL_HANDLER
#define SIGNAL_HANDLER
/* -----------------------------------------------------  
   sighandler.h                                           
   handels to following signals:                          
   SIGIOT         Input/Output problems                   
   SIGFPE         Floating point exceptions               
   SIGBUS         Bus error                               
   SIGSEGV        Segmentation fault                      
   SIGXCPU        CPU time limit exceeded                 
   SIGXFSZ        File size limit exceeded                
   if most of those exception are encountered the system  
   tries to exit gracefully, but with some it dies        
   anyway, but tries to say it why in a way which is     
   for humans better understandable                       
   -----------------------------------------------------  
   part of the lamarc package                             

   P. Beerli                                              
   ----------------------------------------------------- */
#ifdef NEXTAPP
#include <libc.h>
#endif
#define ON 1
#define OFF 0
extern void signalhandling (long switcher);
extern void signalhandler (int sig);
#define error(x) sig_error(x, __FILE__,__LINE__)
extern void sig_error (char string[], char filename[], long line);
#ifdef NEXTAPP
extern void malloc_error_found (int value);
#endif
//#ifdef LAMARC_MALLOC
//extern void *LAMARC_malloc (const long size, const char file[], const long line);
//extern void *LAMARC_calloc (const long repeats, const long size, const char file[], const long line);
//extern void *LAMARC_realloc (void *ptr, const long size, const char file[], const long line);

//#define malloc(a) LAMARC_malloc(a,__FILE__,__LINE__)
//#define calloc(b,a) LAMARC_calloc(b, a,__FILE__,__LINE__)
//#define realloc(c,a) LAMARC_realloc(c, a,__FILE__,__LINE__)
//#endif
#endif
