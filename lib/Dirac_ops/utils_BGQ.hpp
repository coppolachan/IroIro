/*!--------------------------------------------------------------------------
 * @file utils_BGQ.hpp
 * @brief typedef & grobal functions to be used on BGQ
 *Time-stamp: <2013-07-07 08:43:34 noaki>
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
typedef struct GaugeConf{
  double _Complex v[9];
}GaugePtr;

#define ENABLE_THREADING

double omp_norm(Spinor* pointer,int is,int ns,int nid,int tid){
  double tSum;
  BGWilsonLA_Norm(&tSum,pointer+is,ns);
  tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
  if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
  tSum = BGQThread_ScatterDouble(tSum,0,tid,nid);
  return sqrt(tSum);
}

#endif
