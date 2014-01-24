/*!--------------------------------------------------------------------------
 * @file utils_BGQ.cpp

 * Time-stamp: <2013-08-23 11:42:43 cossu>
 *-------------------------------------------------------------------------*/

#include "utils_BGQ.hpp"

double omp_norm(Spinor* pointer,int is,int ns,int nid,int tid){
  double tSum;
  BGWilsonLA_Norm(&tSum,pointer+is,ns);
  tSum = BGQThread_GatherDouble(tSum,0,tid,nid);
  if(tid == 0) tSum = Communicator::instance()->reduce_sum(tSum);
  tSum = BGQThread_ScatterDouble(tSum,0,tid,nid);
  return sqrt(tSum);
}
