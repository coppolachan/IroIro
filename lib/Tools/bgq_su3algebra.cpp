/*--------------------------------------------------------------------
        QCD matrix kernel library

        Copyright 2009-2012 IBM Research - Tokyo, IBM Corporation

        IBM Confidential
---------------------------------------------------------------------*/
#include "include/bgq_su3algebra.h"
#include <complex.h>

void BGWilsonMatLA_Add(QCDComplex* pV,QCDComplex* pW,int ns)
{
  register double _Complex* pV0;
  register double _Complex* pW0;
        int i,j;

        pV0 = pV;
        pW0 = pW;

        for(i=0;i<ns;i++){
                for(j=0;j<QCDLA_MAT;j++){
                        *(pV0 + j) += *(pW0 + j);
                }
                pV0 += QCDLA_MAT;
                pW0 += QCDLA_MAT;
        }
}

void BGWilsonMatLA_Sub(QCDComplex* pV,QCDComplex* pW,int ns)
{
  register double _Complex* pV0;
  register double _Complex* pW0;
        int i,j;

        pV0 = pV;
        pW0 = pW;

        for(i=0;i<ns;i++){
                for(j=0;j<QCDLA_MAT;j++){
                        *(pV0 + j) -= *(pW0 + j);
                }
                pV0 += QCDLA_MAT;
                pW0 += QCDLA_MAT;
        }
}

void BGWilsonMatLA_Add_ND(QCDComplex* pV,QCDComplex* pW,int ns)
{
  register double _Complex* pV0;
  register double _Complex* pW0;

  pV0 = pV;
  pW0 = pW;
  
  for(int i=0; i<ns; ++i){
    *(pV0 +0) += conj(*(pW0 +0));
    *(pV0 +1) += conj(*(pW0 +3));
    *(pV0 +2) += conj(*(pW0 +6));
    *(pV0 +3) += conj(*(pW0 +1));
    *(pV0 +4) += conj(*(pW0 +4));
    *(pV0 +5) += conj(*(pW0 +7));
    *(pV0 +6) += conj(*(pW0 +2));
    *(pV0 +7) += conj(*(pW0 +5));
    *(pV0 +8) += conj(*(pW0 +8));
    
    pV0 += QCDLA_MAT;
    pW0 += QCDLA_MAT;
  }
}


void BGWilsonLA_MatMultScalar(QCDComplex* pV,double PRF,int ns)
{
        register double _Complex* pV0;
        int i,j;

        pV0 = pV;

        for(i=0;i<ns;i++){
                for(j=0;j<QCDLA_MAT;j++){
                        *(pV0 + j) *= PRF;
                }
                pV0 += QCDLA_MAT;
        }
}

void BGWilsonLA_MatEquate(QCDComplex* pV,QCDComplex* pW,int ns)
{
        register double _Complex* pV0;
        register double _Complex* pW0;
        int i,j;

        pV0 = pV;
        pW0 = pW;

        for(i=0;i<ns;i++){
                for(j=0;j<QCDLA_MAT;j++){
                        *(pV0 + j) = *(pW0 + j);
                }
                pV0 += QCDLA_MAT;
                pW0 += QCDLA_MAT;
        }
}


void BGWilsonLA_MatZero(QCDComplex* pV,int ns)
{
        register double _Complex* pV0;
	int i,j;

        pV0 = pV;
        
        for(i=0;i<ns;i++){
                for(j=0;j<QCDLA_MAT;j++){
                        *(pV0 + j) = 0;
                }
                pV0 += QCDLA_MAT;
        }
}

void BGWilsonLA_MatUnity(QCDComplex* pV,int ns)
{
        register double _Complex* pV0;
	int i,j;

        pV0 = pV;
        
        for(i=0;i<ns;i++){
	  for(j=0;j<QCD_COLOR;j++){
	    *(pV0+j*QCD_COLOR+j) = 1.0;
	  }
	  pV0 += QCDLA_MAT;
        }
}



