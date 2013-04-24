//-----------------------------------------------------------------------
// communicator_dummy.cpp
//-----------------------------------------------------------------------
#include "communicator.hpp"
#include "comm_io.hpp"
#include <iostream>

using namespace std;

int Communicator::my_rank_;
int Communicator::Nproc_;
int Communicator::ipe_[Ndim_max_];
int Communicator::nd_up_[Ndim_max_];
int Communicator::nd_dn_[Ndim_max_];

Communicator* Communicator::instance(){
  static Communicator instance_;
  return &instance_;
}

void Communicator::setup(){
  CCIO::cout.init(&std::cout);
  CCIO::cerr.init(&std::cout);
  my_rank_=0;
  Nproc_=1;

  cout << "Communicator using single node mode\n";
}

//----------------------------------------------------------
Communicator::~Communicator(){}
//----------------------------------------------------------


int Communicator::nodeid(int,int,int,int)const{ return my_rank_;}

bool Communicator::primaryNode(){return true;}

void Communicator::transfer_fw(double *bin,double *data,int size,int)const{
  for(int i=0; i<size; ++i) bin[i] = data[i];
}
void Communicator::transfer_fw(valarray<double>& bin, 
			       const valarray<double>& data,int)const{
  bin = data;
}
void Communicator::transfer_fw(valarray<double>& bin, 
			       const valarray<double>& data,
			       const vector<int>& index,int)const{
  for(int i=0; i<index.size(); ++i) bin[i] = data[index[i]];
}

void Communicator::transfer_bk(double *bin,double *data,int size,int)const{
  for(int i=0; i<size; ++i) bin[i] = data[i];
}
void Communicator::transfer_bk(valarray<double>& bin, 
			       const valarray<double>& data,int)const{
  bin = data;
}
void Communicator::transfer_bk(valarray<double>& bin, 
			       const valarray<double>& data,
			       const vector<int>& index,int)const{
  for(int i=0; i<index.size(); ++i) bin[i] = data[index[i]];
}

void Communicator::send_1to1(double *bin,
			     double *data,int size,int,int,int)const{
  for(int i=0; i<size; ++i) bin[i] = data[i];
}

void Communicator::send_1to1(valarray<double>& bin, 
			     const valarray<double>& data,int,int,int,int)const{
  bin = data;
}

double Communicator::reduce_sum(double a)const{ return a;}
void Communicator::sync()const{}
void Communicator::broadcast(int size, int &data, int sender)const{}
void Communicator::broadcast(int size, double &data, int sender)const{}
void Communicator::allgether(double *bin,double *data,int size)const{
  bin = data;}
void Communicator::allgether(valarray<double>& bin,
			     const valarray<double>& data)const{ bin = data;}
int Communicator::reduce_max(double& val,int& idx,int size)const{return 0;}
int Communicator::reduce_min(double& val,int& idx,int size)const{return 0;}
