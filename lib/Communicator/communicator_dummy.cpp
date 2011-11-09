//-----------------------------------------------------------------------
// communicator_dummy.cpp
//-----------------------------------------------------------------------
#include "Communicator/communicator.h"
#include "include/commonPrms.h"
#include "include/field.h"
#include "comm_io.hpp"
#include <stdio.h>
#include <iostream>
#include <cstdarg>
#include <cassert>
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
  int nodeid = my_rank_;

  int NPEx = CommonPrms::instance()->NPEx();
  int NPEy = CommonPrms::instance()->NPEy();
  int NPEz = CommonPrms::instance()->NPEz();
  int NPEt = CommonPrms::instance()->NPEt();

  int px =  nodeid %NPEx;
  int py = (nodeid/NPEx) %NPEy;
  int pz = (nodeid/(NPEx*NPEy)) %NPEz;
  int pt = (nodeid/(NPEx*NPEy*NPEz)) %NPEt;

  ipe_[0] = px;
  ipe_[1] = py;
  ipe_[2] = pz;
  ipe_[3] = pt;

  nd_up_[0] = ((px+1)%NPEx) +py*NPEx +pz*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_up_[1] = px +((py+1)%NPEy)*NPEx +pz*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_up_[2] = px +py*NPEx +((pz+1)%NPEz)*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_up_[3] = px +py*NPEx +pz*NPEx*NPEy +((pt+1)%NPEt)*NPEx*NPEy*NPEz;

  nd_dn_[0] = ((px-1+NPEx)%NPEx) +py*NPEx +pz*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_dn_[1] = px +((py-1+NPEy)%NPEy)*NPEx +pz*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_dn_[2] = px +py*NPEx +((pz-1+NPEz)%NPEz)*NPEx*NPEy +pt*NPEx*NPEy*NPEz;
  nd_dn_[3] = px +py*NPEx +pz*NPEx*NPEy +((pt-1+NPEt)%NPEt)*NPEx*NPEy*NPEz;

  CCIO::cout.init(&std::cout);
}

//----------------------------------------------------------
Communicator::~Communicator(){}
//----------------------------------------------------------
bool Communicator::primaryNode(){return true;}

int Communicator::nodeid(int x,int y,int z,int t){
  int NPEx = CommonPrms::instance()->NPEx();
  int NPEy = CommonPrms::instance()->NPEy();
  int NPEz = CommonPrms::instance()->NPEz();
  int NPEt = CommonPrms::instance()->NPEt();
  return x+NPEx*(y+NPEy*(z+NPEz*t));
}
//-------------------------------------------------------------------
void Communicator::transfer_fw(double *bin, double *data,
			       int size,int){
  for(int i=0; i<size; i++) bin[i] = data[i];
}
//--------------------------------------------------------------------
void Communicator::transfer_bk(double *bin, double *data,
			       int size,int){
  for(int i=0; i<size; i++) bin[i] = data[i];
}
//--------------------------------------------------------------------
void Communicator::transfer_fw(valarray<double>& bin, 
			       const valarray<double>& data, 
			       int size,int dir){
  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];

  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_fw(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin[i] = abin[i];

  delete[] abin;
  delete[] adata;
}

//--------------------------------------------------------------------
void Communicator::transfer_bk(valarray<double>& bin, 
			       const valarray<double>& data,
			       int size,int dir){
  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_bk(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin[i] = abin[i];

  delete[] abin;
  delete[] adata;

}

//---------------------------------------------------------------------
void Communicator::transfer_fw(Field& bin, const Field& data,
			       int size, int dir){

  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];

  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_fw(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin.set(i, abin[i]);

  delete[] abin;
  delete[] adata;

}
//---------------------------------------------------------------------
void Communicator::transfer_bk(Field& bin, const Field& data, 
			       int size, int dir){
  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];

  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_bk(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin.set(i, abin[i]);

  delete[] abin;
  delete[] adata;

}
//---------------------------------------------------------------------
void Communicator::send_1to1(double *bin,  double *data, int size,
			     int p_to, int p_from, int tag){
  for(int i=0; i<size; ++i) bin[i] = data[i];
}

//---------------------------------------------------------------------
void Communicator::send_1to1(valarray<double>& bin, 
			     const valarray<double>& data, int size,
			     int p_to, int p_from, int tag){
  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];

  for(int i=0; i<size; ++i) adata[i] = data[i];
  send_1to1(abin,adata,size,p_to,p_from,tag);
  for(int i=0; i<size; ++i) bin[i] = abin[i];

  delete[] abin;
  delete[] adata;

}
//---------------------------------------------------------------------
void Communicator::send_1to1(Field& bin, const Field& data, int size,
			     int p_to, int p_from, int tag){

  double *abin, *adata;

  abin  = new double[size];
  adata = new double[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  send_1to1(abin,adata,size,p_to,p_from,tag);
  for(int i=0; i<size; ++i) bin.set(i, abin[i]);

  delete[] abin;
  delete[] adata;

}
//--------------------------------------------------------------------
double Communicator::reduce_sum(double& a){ return a;}

//--------------------------------------------------------------------
void Communicator::sync(){
  //    MPI_Barrier(MPI_COMM_WORLD);
}

//--------------------------------------------------------------------
void Communicator::broadcast(int size, int &data, int sender){
  //  MPI_Bcast( &data, size, MPI_INT, sender, MPI_COMM_WORLD );
}

//--------------------------------------------------------------------
void Communicator::broadcast(int size, double &data, int sender){
  //  MPI_Bcast( &data, size, MPI_DOUBLE, sender, MPI_COMM_WORLD );
}

//**************************************************END*****

