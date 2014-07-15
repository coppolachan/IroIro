/*!
 * @file communicator_single.cpp
 *
 * @brief Definition of single node environment Communicator classes
 *
 * Time-stamp: <2014-07-15 11:57:59 neo>
 *
 */

#include "communicator.hpp"
#include "comm_io.hpp"

// Necessary for the static variables
int Communicator::my_rank_;
int Communicator::Nproc_;
int Communicator::ipe_[Ndim_max_];
int Communicator::nd_up_[Ndim_max_];
int Communicator::nd_dn_[Ndim_max_];

//----------------------------------------------------------
Communicator::~Communicator(){}
//----------------------------------------------------------

Communicator* Communicator::instance(){
  static Communicator instance_;
  return &instance_;
}

void Communicator::setup(){
  CCIO::cout.init(&std::cout);
  CCIO::cerr.init(&std::cout);

  // Trivial initializations
  my_rank_ = 0;
  Nproc_   = 1;

  std::cout << "Communicator using single node mode\n";
}


int Communicator::nodeid(int,int,int,int)const{ return 0;}

bool Communicator::primaryNode(){return true;}

void Communicator::transfer_fw(double *bin,
			       double *data,
			       int size,int)const{
  for(int i=0;  i < size; ++i) bin[i] = data[i];
}
void Communicator::transfer_fw(varray_double& bin, 
			       const varray_double& data,int)const{
  bin = data;
}
void Communicator::transfer_fw(varray_double& bin, 
			       const varray_double& data,
			       const vector_int& index,int)const{
  for(int i=0; i<index.size(); ++i) bin[i] = data[index[i]];
}

void Communicator::transfer_bk(double *bin,double *data,int size,int)const{
  for(int i=0; i<size; ++i) bin[i] = data[i];
}
void Communicator::transfer_bk(varray_double& bin, 
			       const varray_double& data,int)const{
  bin = data;
}
void Communicator::transfer_bk(varray_double& bin, 
			       const varray_double& data,
			       const vector_int& index,int)const{
  for(int i=0; i<index.size(); ++i) bin[i] = data[index[i]];
}

void Communicator::send_1to1(double *bin,
			     double *data,int size,
			     int,int,int)const{
  for(int i=0; i<size; ++i) bin[i] = data[i];
}

void Communicator::send_1to1(varray_double& bin, 
			     const varray_double& data,
			     int,int,int,int)const{
  bin = data;
}

void Communicator::sync()const{}
void Communicator::broadcast(int size, 
			     int &data, 
			     int sender)const{}

void Communicator::broadcast(int size, 
			     double &data, 
			     int sender)const{}

void Communicator::allgather(double *bin,
			     double *data,
			     int size)const{
  bin = data;
}

void Communicator::allgather(varray_double& bin,
			     const varray_double& data)const{ bin = data;}

int Communicator::reduce_max(double& val,int& idx,int size)const{return 0;}

int Communicator::reduce_min(double& val,int& idx,int size)const{return 0;}

double Communicator::reduce_sum(double a)const{ return a;}

uint32_t Communicator::reduce_sum(uint32_t a)const{ return a;}
uint64_t Communicator::reduce_sum(uint64_t a)const{ return a;}
