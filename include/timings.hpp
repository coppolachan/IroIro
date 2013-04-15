/*
  @file: timing.hpp

  Timing information routines for BGQ


 */


#ifndef TIMING_BGQ_HPP_
#define TIMING_BGQ_HPP_

//now empty
#define FINE_TIMING_START(t)
#define FINE_TIMING_END(t)

#ifdef IBM_BGQ_WILSON

#define BGQ_HERTZ 1600000000.0

#include <spi/include/kernel/process.h> //for GetTimeBase


namespace timingBGQ {

  long double Timer(long double t = 0) {
    long double new_t = ((long double)GetTimeBase())/BGQ_HERTZ; 
    
    if (t == 0)
      return new_t;
    else 
      return new_t - t;
  }

  // Redefine the timing routines in case of BGQ 
#undef FINE_TIMING_START
#define FINE_TIMING_START(t) t = 0; \ 
  t = timingBGQ::Timer()

#undef  FINE_TIMING_END
#define FINE_TIMING_END(t) t = timingBGQ::Timer(t)

}

#endif

#endif
