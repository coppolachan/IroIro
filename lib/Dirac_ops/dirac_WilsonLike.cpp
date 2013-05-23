/*! @file dirac_WilsonLike.cpp
 *  @brief implementation of member functions of Dirac_WilsonLike
 Time-stamp: <2013-05-22 21:12:16 noaki>
 */

#include "dirac_WilsonLike.hpp"

void DiracWilsonLike::gammaXcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    w[r0(c)] = f[i3(c)];  w[i0(c)] =-f[r3(c)];
    w[r1(c)] = f[i2(c)];  w[i1(c)] =-f[r2(c)];
    w[r2(c)] =-f[i1(c)];  w[i2(c)] = f[r1(c)];
    w[r3(c)] =-f[i0(c)];  w[i3(c)] = f[r0(c)];
  }
}
void DiracWilsonLike::gammaYcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    w[r0(c)] =-f[r3(c)];  w[i0(c)] =-f[i3(c)];
    w[r1(c)] = f[r2(c)];  w[i1(c)] = f[i2(c)];
    w[r2(c)] = f[r1(c)];  w[i2(c)] = f[i1(c)];
    w[r3(c)] =-f[r0(c)];  w[i3(c)] =-f[i0(c)];
  }
}
void DiracWilsonLike::gammaZcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    w[r0(c)] = f[i2(c)];  w[i0(c)] =-f[r2(c)];
    w[r1(c)] =-f[i3(c)];  w[i1(c)] = f[r3(c)];
    w[r2(c)] =-f[i0(c)];  w[i2(c)] = f[r0(c)];
    w[r3(c)] = f[i1(c)];  w[i3(c)] =-f[r1(c)];
  }
}
void DiracWilsonLike::gammaTcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    w[r0(c)] = f[r0(c)];  w[i0(c)] = f[i0(c)];
    w[r1(c)] = f[r1(c)];  w[i1(c)] = f[i1(c)];
    w[r2(c)] =-f[r2(c)];  w[i2(c)] =-f[i2(c)];
    w[r3(c)] =-f[r3(c)];  w[i3(c)] =-f[i3(c)];
  }
}
void DiracWilsonLike::gamma5core(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    w[r0(c)] = f[r2(c)];  w[i0(c)] = f[i2(c)];
    w[r1(c)] = f[r3(c)];  w[i1(c)] = f[i3(c)];
    w[r2(c)] = f[r0(c)];  w[i2(c)] = f[i0(c)];
    w[r3(c)] = f[r1(c)];  w[i3(c)] = f[i1(c)];
  }
}

void DiracWilsonLike::projPcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    double fup_r = 0.5*(f[r0(c)] +f[r2(c)]);
    double fup_i = 0.5*(f[i0(c)] +f[i2(c)]);
    double fdn_r = 0.5*(f[r1(c)] +f[r3(c)]);
    double fdn_i = 0.5*(f[i1(c)] +f[i3(c)]);
    w[r0(c)] = fup_r;   w[i0(c)] = fup_i;
    w[r1(c)] = fdn_r;   w[i1(c)] = fdn_i;
    w[r2(c)] = fup_r;   w[i2(c)] = fup_i;
    w[r3(c)] = fdn_r;   w[i3(c)] = fdn_i;
  }
}

void DiracWilsonLike::projMcore(double* w,const double* f)const{
  for(int c=0; c<NC_; ++c){
    double fup_r = 0.5*(f[r0(c)] -f[r2(c)]);
    double fup_i = 0.5*(f[i0(c)] -f[i2(c)]);
    double fdn_r = 0.5*(f[r1(c)] -f[r3(c)]);
    double fdn_i = 0.5*(f[i1(c)] -f[i3(c)]);
    w[r0(c)] = fup_r;   w[i0(c)] = fup_i;
    w[r1(c)] = fdn_r;   w[i1(c)] = fdn_i;
    w[r2(c)] =-fup_r;   w[i2(c)] =-fup_i;
    w[r3(c)] =-fdn_r;   w[i3(c)] =-fdn_i;
  }
}

void Dirac_optimalDomainWall_4D::
BprojCore(double* f,const double* f1,const double* fN)const{ 
  for(int c=0; c<NC_; ++c){
    double fupNr = 0.5*(fN[r0(c)] +fN[r2(c)]);
    double fupNi = 0.5*(fN[i0(c)] +fN[i2(c)]);
    double fdnNr = 0.5*(fN[r1(c)] +fN[r3(c)]);
    double fdnNi = 0.5*(fN[i1(c)] +fN[i3(c)]);

    double fup1r = 0.5*(f1[r0(c)] -f1[r2(c)]);
    double fup1i = 0.5*(f1[i0(c)] -f1[i2(c)]);
    double fdn1r = 0.5*(f1[r1(c)] -f1[r3(c)]);
    double fdn1i = 0.5*(f1[i1(c)] -f1[i3(c)]);

    f[r0(c)] = fupNr +fup1r;   f[i0(c)] = fupNi +fup1i;
    f[r1(c)] = fdnNr +fdn1r;   f[i1(c)] = fdnNi +fdn1i;
    f[r2(c)] = fupNr -fup1r;   f[i2(c)] = fupNi -fup1i;
    f[r3(c)] = fdnNr -fdn1r;   f[i3(c)] = fdnNi -fdn1i;
  }
}

void Dirac_optimalDomainWall_4D::
BprojCore_dag(double* f1,double* fN,const double* f) const{

  // f1 = f5(0), fN = f5(N5_-1)
  for(int c=0; c<NC_; ++c){
    double fup_r = 0.5*(f[r0(c)] +f[r2(c)]);
    double fup_i = 0.5*(f[i0(c)] +f[i2(c)]);
    double fdn_r = 0.5*(f[r1(c)] +f[r3(c)]);
    double fdn_i = 0.5*(f[i1(c)] +f[i3(c)]);

    fN[r0(c)] = fup_r;   fN[i0(c)] = fup_i;
    fN[r1(c)] = fdn_r;   fN[i1(c)] = fdn_i;
    fN[r2(c)] = fup_r;   fN[i2(c)] = fup_i;
    fN[r3(c)] = fdn_r;   fN[i3(c)] = fdn_i;

    fup_r -= f[r2(c)]; //0.5*(f[r0(c)] -f[r2(c)])
    fup_i -= f[i2(c)]; //0.5*(f[i0(c)] -f[i2(c)])
    fdn_r -= f[r3(c)]; //0.5*(f[r1(c)] -f[r3(c)])
    fdn_i -= f[i3(c)]; //0.5*(f[i1(c)] -f[i3(c)])

    f1[r0(c)] = fup_r;   f1[i0(c)] = fup_i;
    f1[r1(c)] = fdn_r;   f1[i1(c)] = fdn_i;
    f1[r2(c)] =-fup_r;   f1[i2(c)] =-fup_i;
    f1[r3(c)] =-fdn_r;   f1[i3(c)] =-fdn_i;
  }
}




