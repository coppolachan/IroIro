/*! 
  @file commonPrms.cpp  
  @brief Defines CommonPrms class
*/

#include "include/commonPrms.h"
#include "include/messages_macros.hpp"
#include<iostream>
#include<cstdlib>

using namespace std;

int CommonPrms::Lx_;
int CommonPrms::Ly_;
int CommonPrms::Lz_;
int CommonPrms::Lt_;
int CommonPrms::Lvol_;

int CommonPrms::Nx_;
int CommonPrms::Ny_;
int CommonPrms::Nz_;
int CommonPrms::Nt_;
int CommonPrms::Nvol_;

int CommonPrms::NPEx_;
int CommonPrms::NPEy_;
int CommonPrms::NPEz_;
int CommonPrms::NPEt_;
int CommonPrms::NP_;

vector<int> CommonPrms::Lsize_(NDIM_);
vector<int> CommonPrms::Nsize_(NDIM_); 

CommonPrms* CommonPrms::instance_= NULL;

CommonPrms* CommonPrms::instance(const Lattice& latt){
  if(instance_== NULL) instance_= new CommonPrms(latt);
  return instance_;
}

CommonPrms* CommonPrms::instance(){
  if(instance_== NULL){
    exit(EXIT_FAILURE);
  }
  return instance_;
}

CommonPrms::CommonPrms(const Lattice& latt){
  if(latt.stdinput){
    cin >> Lx_>> Ly_>> Lz_>> Lt_;
    NPEx_ = 1;
    NPEy_ = 1;
    NPEz_ = 1;
    NPEt_ = 1;
  }else{
    Lx_= latt.Lx;
    Ly_= latt.Ly;
    Lz_= latt.Lz;
    Lt_= latt.Lt;
    
    NPEx_= latt.NPEx;
    NPEy_= latt.NPEy;
    NPEz_= latt.NPEz;
    NPEt_= latt.NPEt;
  }

  Lvol_= Lx_*Ly_*Lz_*Lt_;
  NP_= NPEx_*NPEy_*NPEz_*NPEt_;

  Nx_= Lx_/NPEx_;
  Ny_= Ly_/NPEy_;
  Nz_= Lz_/NPEz_;
  Nt_= Lt_/NPEt_;
  Nvol_= Lvol_/NP_;

  Lsize_.push_back(Lx_);
  Lsize_.push_back(Ly_);
  Lsize_.push_back(Lz_);
  Lsize_.push_back(Lt_);

  Nsize_.push_back(Nx_);
  Nsize_.push_back(Ny_);
  Nsize_.push_back(Nz_);
  Nsize_.push_back(Nt_);
}
