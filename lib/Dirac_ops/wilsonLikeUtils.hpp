#ifndef WILSONLIKEUTILS_INCLUDED
#define WILSONLIKEUTILS_INCLUDED

#include "include/macros.hpp"

class GammaMatrix{
  int Ncol_;
  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(Ncol_+c);}
  int r2(int c)const{return 2*(2*Ncol_+c);}
  int r3(int c)const{return 2*(3*Ncol_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(Ncol_+c)+1;}
  int i2(int c)const{return 2*(2*Ncol_+c)+1;}
  int i3(int c)const{return 2*(3*Ncol_+c)+1;}

public:
  GammaMatrix(int Ncol=NC_):Ncol_(Ncol){}

  void gammaXcore(double*,const double*)const;
  void gammaYcore(double*,const double*)const;
  void gammaZcore(double*,const double*)const;
  void gammaTcore(double*,const double*)const;
  void gamma5core(double*,const double*)const;
  void projPcore(double*,const double*)const;
  void projMcore(double*,const double*)const;
  
  static void (GammaMatrix::*gamma[])(double*,const double*)const;
};

class DW5dMatrix{
  int Ncol_;
  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(Ncol_+c);}
  int r2(int c)const{return 2*(2*Ncol_+c);}
  int r3(int c)const{return 2*(3*Ncol_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(Ncol_+c)+1;}
  int i2(int c)const{return 2*(2*Ncol_+c)+1;}
  int i3(int c)const{return 2*(3*Ncol_+c)+1;}
public:
  DW5dMatrix(int Ncol=NC_):Ncol_(Ncol){}
  
  void BprojCore(double* f1,const double* fN,const double* f) const;
  void BprojCore_dag(double* f1,double* fN,const double* f) const;
};

#endif
