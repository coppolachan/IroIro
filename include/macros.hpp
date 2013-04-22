/*!
 * @file macros.hpp
 *
 * @brief Include several useful macros and definitions
 *
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_

#include <sys/time.h>

#ifdef HAVE_CONFIG_H
#include "iroiro_config.h"
#endif

#include <iostream>
#include <sstream>
#include <vector>

// Predefined verbosity levels
#define BASE_VERB_LEVEL 1
#define ACTION_VERB_LEVEL 1
#define TIMING_VERB_LEVEL 2
#define SOLV_ITER_VERB_LEVEL 3
#define SOLV_MONITOR_VERB_LEVEL 1
#define DEBUG_VERB_LEVEL 5

#define NC_ 3  //Number of colours, now fixed here
#define ND_ 4  //Number of dirac spinors, now fixed here
#define NDIM_ 4  //Number of lattice dimensions, now fixed here

#define TIMING_START				\
  timeval start, end;				\
  gettimeofday(&start,NULL)

#define TIMING_END(var)							\
  gettimeofday(&end,NULL);						\
  var = (end.tv_sec - start.tv_sec) * 1000.0;				\
  var = var + (end.tv_usec - start.tv_usec) / 1000.0   // us to ms

class NullType{};


typedef std::vector<int> vector_int;
typedef std::vector<double> vector_double;



#endif
