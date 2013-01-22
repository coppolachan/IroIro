/*--------------------------------------------------------------------
        QCD matrix kernel library

        Copyright 2009-2012 IBM Research - Tokyo, IBM Corporation

        IBM Confidential
---------------------------------------------------------------------*/

#ifndef BGQ_SU3ALGEBRA_H_
#define BGQ_SU3ALGEBRA_H_

#ifdef IBM_BGQ_WILSON
#include "/bghome/scbadm/ibm-doi/lib/lib_wilson/qcd.h"

#define QCDLA_MAT          9

void BGWilsonMatLA_Add(QCDComplex* pV,QCDComplex* pW,int ns);
void BGWilsonLA_MatMultScalar(QCDComplex* pV,double PRF,int ns);
void BGWilsonLA_MatEquate(QCDComplex* pV,QCDComplex* pW,int ns);

#endif

#endif
