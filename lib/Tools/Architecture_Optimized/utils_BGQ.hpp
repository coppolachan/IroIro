/*!--------------------------------------------------------------------------
 * @file utils_BGQ.hpp
 * @brief typedef & global functions to be used on BGQ

 * Time-stamp: <2014-04-17 12:49:49 neo>
 *-------------------------------------------------------------------------*/
#ifndef UTILS_BGQ_INCLUDED
#define UTILS_BGQ_INCLUDED

#include "include/timings.hpp"
#include "bgqwilson.h"
#include "bgqthread.h"
#include <omp.h>
#include <hwi/include/bqc/A2_inlines.h>
#include "include/messages_macros.hpp"

typedef struct FermionSpinor{
  double _Complex v[12];
}Spinor;
typedef struct GaugeConfig{
  double _Complex v[9];
}GaugePtr;

#define ENABLE_THREADING

double omp_norm(Spinor* pointer,int is,int ns,int nid,int tid);

void MonitorBGQMemory();


#endif
