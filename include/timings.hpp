/*
  @file: timing.hpp

  Timing information routines for BGQ


 */


#ifndef TIMING_BGQ_HPP_
#define TIMING_BGQ_HPP_


//now empty
#define FINE_TIMING_START(t)
#define FINE_TIMING_END(t)

#if (VERBOSITY>=TIMING_VERB_LEVEL)
#ifdef IBM_BGQ_WILSON
// Redefine the timing routines in case of BGQ 
// and only if the fime timing has been enabled
// VERBOSITY>=TIMING_VERB_LEVEL

namespace timingBGQ {
  long double Timer(long double t = 0);
}

/////////////////////////////////////////////
#undef FINE_TIMING_START
#define FINE_TIMING_START(t) t = 0; \
                             t = timingBGQ::Timer()
  
#undef  FINE_TIMING_END
#define FINE_TIMING_END(t) t = timingBGQ::Timer(t)
///////////////////////////////////////////
  
#endif // IBM_BGQ_WILSON

#endif  // (VERBOSITY >= TIMING_VERB_LEVEL)

#endif // TIMING_BGQ_HPP_
