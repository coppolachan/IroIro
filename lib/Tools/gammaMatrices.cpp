#include "gammaMatrices.hpp"

namespace GammaMatrices{
  
  GammaResult Gamma::operator() (int s) const{
    GammaResult res;
    res.spn = idx_[s];
    res.facr = fac_[2*s];
    res.faci = fac_[2*s+1];
    return res;
  }
  
  Unit::Unit(){
    int idx[] = {0,1,2,3,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0;
    fac_[2] = 1.0;
    fac_[4] = 1.0;
    fac_[6] = 1.0;
  }

  Gamma1::Gamma1(){
    int idx[] = {3,2,1,0,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] =-1.0;
    fac_[3] =-1.0;
    fac_[5] = 1.0;
    fac_[7] = 1.0;
  }

  Gamma2::Gamma2(){
    int idx[] = {3,2,1,0,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] =-1.0; 
    fac_[2] = 1.0; 
    fac_[4] = 1.0; 
    fac_[6] =-1.0; 
  }

  Gamma3::Gamma3(){
    int idx[] = {2,3,0,1,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] =-1.0;
    fac_[3] = 1.0;
    fac_[5] = 1.0;
    fac_[7] =-1.0;
  }

  Gamma4::Gamma4(){
    int idx[] = {0,1,2,3,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0; 
    fac_[2] = 1.0; 
    fac_[4] =-1.0; 
    fac_[6] =-1.0; 
  }

  Gamma5::Gamma5(){
    int idx[] = {2,3,0,1,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0; 
    fac_[2] = 1.0; 
    fac_[4] = 1.0; 
    fac_[6] = 1.0; 
  }

  Gamma1_5::Gamma1_5(){
    int idx[] = {1,0,3,2,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] =-1.0;
    fac_[3] =-1.0;
    fac_[5] = 1.0;
    fac_[7] = 1.0;
  }

  Gamma2_5::Gamma2_5(){
    int idx[] = {1,0,3,2,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] =-1.0; 
    fac_[2] = 1.0; 
    fac_[4] = 1.0; 
    fac_[6] =-1.0; 
  }

  Gamma3_5::Gamma3_5(){
    int idx[] = {0,1,2,3,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] =-1.0;
    fac_[3] = 1.0;
    fac_[5] = 1.0;
    fac_[7] =-1.0;
  }

  Gamma4_5::Gamma4_5(){
    int idx[] = {2,3,0,1,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0;  
    fac_[2] = 1.0;  
    fac_[4] =-1.0;  
    fac_[6] =-1.0;  
  }

  // [,] denotes the commutation relation.
  Sigma12_5::Sigma12_5(){ // (1/2)*[Gamma1,Gamma2]*Gamma5
    int idx[] = {2,3,0,1,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] = 1.0;
    fac_[3] =-1.0;
    fac_[5] = 1.0;
    fac_[7] =-1.0;
  }

  Sigma13_5::Sigma13_5(){ // (1/2)*[Gamma1,Gamma3]*Gamma5             
    int idx[] = {3,2,1,0,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0;
    fac_[2] =-1.0;
    fac_[4] = 1.0;
    fac_[6] =-1.0;
  }

  Sigma14_5::Sigma14_5(){ // (1/2)*[Gamma1,Gamma4]*Gamma5       
    int idx[] = {1,0,2,3,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0;
    fac_[2] = 1.0;
    fac_[4] = 1.0;
    fac_[6] = 1.0;
  }

  Sigma23_5::Sigma23_5(){ // (1/2)*[Gamma2,Gamma3]*Gamma5
    int idx[] = {3,2,1,0,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] = 1.0;
    fac_[3] = 1.0;
    fac_[5] = 1.0;
    fac_[7] = 1.0;
  }

  Sigma24_5::Sigma24_5(){ // (1/2)*[Gamma2,Gamma4]*Gamma5
    int idx[] = {1,0,3,2,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0;
    fac_[2] = -1.0;
    fac_[4] = 1.0;
    fac_[6] = -1.0;
  }

  Sigma34_5::Sigma34_5(){ // (1/2)*[Gamma3,Gamma4]*Gamma5          
    int idx[] = {0,1,2,3,};
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] = 1.0;
    fac_[3] = -1.0;
    fac_[5] = 1.0;
    fac_[7] = -1.0;
  }

  CConj::CConj(){ //Gamma2*Gamma4
    int idx[] = {3,2,1,0,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0; 
    fac_[2] =-1.0; 
    fac_[4] = 1.0; 
    fac_[6] =-1.0; 
  }

  CGamma1::CGamma1(){ //Gamma2*Gamma4*Gamma1
    int idx[] = {0,1,2,3,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] = 1.0; 
    fac_[3] =-1.0; 
    fac_[5] =-1.0; 
    fac_[7] = 1.0; 
  }

  CGamma2::CGamma2(){ //Gamma2*Gamma4*Gamma2
    int idx[] = {0,1,2,3,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] =-1.0; 
    fac_[2] =-1.0; 
    fac_[4] = 1.0; 
    fac_[6] = 1.0; 
  }

  CGamma3::CGamma3(){ //Gamma2*Gamma4*Gamma3
    int idx[] = {1,0,3,2,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[1] =-1.0; 
    fac_[3] =-1.0; 
    fac_[5] = 1.0; 
    fac_[7] = 1.0; 
  }

  CGamma5::CGamma5(){ //Gamma2*Gamma4*Gamma5
    int idx[] = {1,0,3,2,}; 
    for(int i=0; i<NDIM_; ++i) idx_.push_back(idx[i]);

    fac_[0] = 1.0; 
    fac_[2] =-1.0; 
    fac_[4] = 1.0; 
    fac_[6] =-1.0; 
  }

}
