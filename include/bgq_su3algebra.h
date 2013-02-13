/*--------------------------------------------------------------------
        QCD matrix kernel library

        Copyright 2009-2012 IBM Research - Tokyo, IBM Corporation

        IBM Confidential
---------------------------------------------------------------------*/

#ifndef BGQ_SU3ALGEBRA_H_
#define BGQ_SU3ALGEBRA_H_


#include "/bghome/scbadm/ibm-doi/lib/lib_wilson/qcd.h"

#define QCD_COLOR          3
#define QCDLA_MAT          QCD_COLOR*QCD_COLOR

void BGWilsonMatLA_Add(QCDComplex* pV,QCDComplex* pW,int ns);
void BGWilsonMatLA_Sub(QCDComplex* pV,QCDComplex* pW,int ns);
void BGWilsonMatLA_Add_ND(QCDComplex* pV,QCDComplex* pW,int ns);
void BGWilsonLA_MatMultScalar(QCDComplex* pV,double PRF,int ns);
void BGWilsonLA_MatEquate(QCDComplex* pV,QCDComplex* pW,int ns);
void BGWilsonLA_MatZero(QCDComplex* pV,int ns);
void BGWilsonLA_MatUnity(QCDComplex* pV,int ns);



#endif
