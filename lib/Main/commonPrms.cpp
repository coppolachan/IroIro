/*! 
  @file commonPrms.cpp  
  @brief Defines CommonPrms class
*/

#include "include/commonPrms.h"
#include "include/messages_macros.hpp"
#include<iostream>
#include<cstdlib>

using namespace std;

int CommonPrms::Lvol_;
int CommonPrms::Nvol_;
int CommonPrms::NP_;

vector<int> CommonPrms::Lsize_(NDIM_);
vector<int> CommonPrms::Nsize_(NDIM_); 
vector<int> CommonPrms::Nnodes_(NDIM_); 

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
    cin >> Lsize_[0]>> Lsize_[1]>> Lsize_[2]>> Lsize_[3];

    Nnodes_[0] = 1;
    Nnodes_[1] = 1;
    Nnodes_[2] = 1;
    Nnodes_[3] = 1;
  }else{

    Lsize_[0] = latt.Lx;
    Lsize_[1] = latt.Ly;
    Lsize_[2] = latt.Lz;
    Lsize_[3] = latt.Lt;
    
    Nnodes_[0] = latt.NPEx;
    Nnodes_[1] = latt.NPEy;
    Nnodes_[2] = latt.NPEz;
    Nnodes_[3] = latt.NPEt;
  }
  Lvol_= Lsize_[0]*Lsize_[1]*Lsize_[2]*Lsize_[3];
  NP_= Nnodes_[0]*Nnodes_[1]*Nnodes_[2]*Nnodes_[3];

  Nsize_[0] = Lsize_[0]/Nnodes_[0];
  Nsize_[1] = Lsize_[1]/Nnodes_[1];
  Nsize_[2] = Lsize_[2]/Nnodes_[2];
  Nsize_[3] = Lsize_[3]/Nnodes_[3];

  Nvol_= Lvol_/NP_;
}
