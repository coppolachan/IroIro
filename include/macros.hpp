/*!
 * @file macros.hpp
 *
 * @brief Include several useful macros and definitions
 *
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_

#include <time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// Predefined verbosity levels
#define ACTION_VERB_LEVEL 1
#define SOLV_ITER_VERB_LEVEL 2
#define SOLV_MONITOR_VERB_LEVEL 1
#define DEBUG_VERB_LEVEL 5

#define TIMING_START				\
  clock_t start, end;				\
  start = clock()

#define TIMING_END(var)				\
  end = clock();				\
  var = (end-start)*1000/CLOCKS_PER_SEC



#endif
