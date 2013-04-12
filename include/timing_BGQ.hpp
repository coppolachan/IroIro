/*
  @file: timing.hpp

  Timing information routines for BGQ


 */
//now empty
#define FINE_TIMING_START
#define FINE_TIMING_PRINT(var)

#ifndef TIMING_BGQ_HPP_
#define TIMING_BGQ_HPP_

#ifdef IBM_BGQ_WILSON

#define BGQ_HERTZ 1600000000.0

#include <spi/include/kernel/process.h> //for GetTimeBase


namespace timingBGQ {

  // Returns the Time from the last call
  long double GetTime() {
    static long double old;
    static bool have_time = false;
    long double diff, newest;

    if (have_time) {
      newest = ((long double)GetTimeBase())/BGQ_HERTZ;
      diff = newest -old ;
      old = newest;
      return diff;
    } else {
      old = ((long double)GetTimeBase())/BGQ_HERTZ;
      have_time = true;
      return old;
    }

  }

  // Redefine the timing routines in case of BGQ 
#undef FINE_TIMING_START
#define FINE_TIMING_START timingBGQ::GetTime()

#undef FINE TIMING_END(var)
#define FINE_TIMING_END(var) var = timingBGQ::GetTime()

}

#endif

#endif
