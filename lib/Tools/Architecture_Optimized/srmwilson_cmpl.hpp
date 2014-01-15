#ifndef SRMWILSON_CMPL_INCLUDED
#define SRMWILSON_CMPL_INCLUDED
#include "macros.hpp"

namespace SRCMPL{
  
  inline SRWilsonSU3_MatUnity(double* dat,int size){
    int NCC= 2*NC_*NC_;
    for(int i=0; i<size; ++i){
      for(int cc=0; cc<NCC; ++cc) *(dat +i*NCC +cc) = 0.0;
      for(int dd=0; dd<NC_; ++dd) *(dat +i*NCC +dd*8) = 1.0;
    }
  }

  inline SRWilsonSU3_MatEquate(double* dat1,double* dat2,int size){
    int NCC= 2*NC_*NC_;
    for(int i=0; i<size*NCC; ++i)  *(dat1 +i) = *(dat2 +i);
  }

  inline SRWilsonSU3_MatMultScalar(double* dat,double fac,int size){
    int NCC= 2*NC_*NC_;
    for(int i=0; i<size*NCC; ++i)  *(dat +i) *=  fac;
  }

  inline SRWilsonSU3_MatSub(double* dat,double* sub,int size){
    int NCC= 2*NC_*NC_;
    for(int i=0; i<size*NCC; ++i)  *(dat +i) -= *(sub+i);
  }

  inline SRWilsonSU3_MatAdd(double* dat,double* sub,int size){
    int NCC= 2*NC_*NC_;
    for(int i=0; i<size*NCC; ++i)  *(dat +i) += *(sub+i);
  }
}
#endif
