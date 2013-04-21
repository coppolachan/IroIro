/*
  @file: timing.cpp

  Timing information routines for BGQ


 */

#include "include/timings.hpp"
#include "include/macros.hpp"

#if (VERBOSITY>=TIMING_VERB_LEVEL)
#ifdef IBM_BGQ_WILSON

#define BGQ_HERTZ 1600000000.0
#include <spi/include/kernel/process.h> //for GetTimeBase

namespace timingBGQ{

  long double Timer(long double t = 0) {
    long double new_t = ((long double)GetTimeBase())/BGQ_HERTZ; 
    
    if (t == 0)
      return new_t;
    else 
      return new_t - t;
  }
}


#endif
#endif
