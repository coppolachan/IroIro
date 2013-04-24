/*! @file dirac_WilsonLike.cpp
 *  @brief implementation of member functions of Dirac_WilsonLike
 Time-stamp: <2013-04-23 12:48:37 noaki>
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


