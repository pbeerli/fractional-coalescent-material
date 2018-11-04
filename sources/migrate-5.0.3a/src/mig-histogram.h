#ifndef _MIGHIST_H_
#define _MIGHIST_H_
/*
 * mighistogram material
 *
 
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject
 to the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 * $Id: mig-histogram.h 2157 2013-04-16 22:27:40Z beerli $
 */
#include "migration.h"

#define NOAVERAGE -1.

extern void increase_mighist (mighistloci_fmt * mighistlocus);
extern void setup_mighist (world_fmt * world, option_fmt * options);
extern void print_mighist (world_fmt * world);

extern void minmax (histogram_fmt * hist, float *tempmin, float *tempmax);

extern void destroy_mighist (world_fmt * world);
#endif




