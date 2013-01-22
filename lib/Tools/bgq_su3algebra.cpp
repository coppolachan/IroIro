/*--------------------------------------------------------------------
        QCD matrix kernel library

        Copyright 2009-2012 IBM Research - Tokyo, IBM Corporation

        IBM Confidential
---------------------------------------------------------------------*/

#ifdef IBM_BGQ_WILSON
#include "bgq_su3algebra.h"

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

#endif

