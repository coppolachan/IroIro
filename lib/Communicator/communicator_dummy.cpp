//-----------------------------------------------------------------------
// communicator_dummy.cpp
//-----------------------------------------------------------------------
#include "communicator.h"
//#include <cassert>
using namespace std;

int Communicator::my_rank_;
int Communicator::Nproc_;
int Communicator::ipe_[Ndim_max_];
int Communicator::nd_up_[Ndim_max_];
int Communicator::nd_dn_[Ndim_max_];

Communicator* Communicator::instance(){
  static Communicator* instance_;
  return instance_;
}

void Communicator::setup(){
  my_rank_=0;
  Nproc_=1;
}
int Communicator::nodeid(int,int,int,int)const{ return my_rank_;}

bool Communicator::primaryNode(){return true;}

void Communicator::transfer_fw(double *bin,double *data,int size,int){
  for(int i=0; i<size; i++) bin[i] = data[i];
}
void Communicator::transfer_fw(valarray<double>& bin, 
			       const valarray<double>& data,int,int){
  bin = data;
}

void Communicator::transfer_bk(double *bin,double *data,int size,int){
  for(int i=0; i<size; i++) bin[i] = data[i];
}
void Communicator::transfer_bk(valarray<double>& bin, 
			       const valarray<double>& data,int,int){
  bin = data;
}

void Communicator::send_1to1(double *bin,
			     double *data,int size,int,int,int){
  for(int i=0; i<size; ++i) bin[i] = data[i];
}

void Communicator::send_1to1(valarray<double>& bin, 
			     const valarray<double>& data,int,int,int,int){
  bin = data;
}

double Communicator::reduce_sum(double a){ return a;}
void Communicator::sync(){}
void Communicator::broadcast(int size, int &data, int sender){}
void Communicator::broadcast(int size, double &data, int sender){}


