#ifndef SRMWILSON_CMPL_INCLUDED
#define SRMWILSON_CMPL_INCLUDED
#include "macros.hpp"

namespace SRCMPL{
  
  inline SRWilsonSU3_MatZero(double* dat,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size*CC2; ++i) *(dat+i) = 0.0;
  }
  
  inline SRWilsonSU3_MatUnity(double* dat,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size; ++i){
      for(int cc=0; cc<CC2; ++cc) *(dat +i*CC2 +cc) = 0.0;
      for(int dd=0; dd<NC_; ++dd) *(dat +i*CC2 +dd*8) = 1.0;
    }
  }

  inline SRWilsonSU3_MatEquate(double* dat1,double* dat2,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size*CC2; ++i)  *(dat1 +i) = *(dat2 +i);
  }

  inline SRWilsonSU3_MatMultScalar(double* dat,double fac,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size*CC2; ++i)  *(dat +i) *=  fac;
  }

  inline SRWilsonSU3_MatSub(double* dat,double* sub,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size*CC2; ++i)  *(dat +i) -= *(sub+i);
  }

  inline SRWilsonSU3_MatAdd(double* dat,double* add,int size){
    int CC2= 2*NC_*NC_;
    for(int i=0; i<size*CC2; ++i)  *(dat +i) += *(add+i);
  }

  inline SRWilsonSU3_MatAdd_ND(double* dat,double* add,int size){
    int CC2= 2*NC_*NC_;

    for(int v=0; v<size; ++v){
      for(int a=0; a<NC_; ++a){
	int aa = v*CC2 +2*(NC_*a+a);
	*(dat +aa  ) += *(add +aa  );
	*(dat +aa+1) -= *(add +aa+1);

	for(int b=a+1; b<NC_; ++b){
	  int ab = v*CC2 +2*(NC_*a+b);
	  int ba = v*CC2 +2*(NC_*b+a);

	  *(dat +ab  ) += *(add +ba  );
	  *(dat +ab+1) -= *(add +ba+1);

	  *(dat +ba  ) += *(add +ab  );
	  *(dat +ba+1) -= *(add +ab+1);
	}
      }
    }
  }
}
#endif
