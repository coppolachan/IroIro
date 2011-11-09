//-----------------------------------------------------------------
// communicator_mpi.cpp
//-----------------------------------------------------------------
#include "communicator.h"
#include "mpi.h"
#include "include/commonPrms.h"
#include "include/field.h"
#include "comm_io.hpp"
#include <stdio.h>
#include <iostream>
#include <cstdarg>
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
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &Nproc_);

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
  CCIO::cerr.init(&std::cerr);

  if(my_rank_==0) cout << "Communicator using MPI with "<< Nproc_ << " processes.\n";
}
//----------------------------------------------------------
Communicator::~Communicator(){
  MPI_Finalize();
}

bool Communicator::primaryNode(){
  if (my_rank_==0) 
    return true;
  else
    return false;
}

//----------------------------------------------------------
int Communicator::nodeid(int x, int y, int z, int t){
  int NPEx = CommonPrms::instance()->NPEx();
  int NPEy = CommonPrms::instance()->NPEy();
  int NPEz = CommonPrms::instance()->NPEz();
  int NPEt = CommonPrms::instance()->NPEt();
  return x +NPEx*(y +NPEy*(z +NPEz*t));
}
  
//----------------------------------------------------------
void Communicator::transfer_fw(double *bin, double *data,int size,int dir){
  MPI_Status  status;
  int p_send = nd_dn_[dir];
  int p_recv = nd_up_[dir];

  int tag1 = dir*Nproc_+my_rank_;
  int tag2 = dir*Nproc_+p_recv;

  MPI_Sendrecv( data, size, MPI_DOUBLE, p_send, tag1,
		bin,  size, MPI_DOUBLE, p_recv, tag2,
		MPI_COMM_WORLD, &status );
}
//----------------------------------------------------------
void Communicator::transfer_bk(double *bin,double *data,int size,int dir){
  MPI_Status  status;
  int p_send = nd_up_[dir];
  int p_recv = nd_dn_[dir];
  int Ndim = CommonPrms::instance()->Ndim();
  int tag1 = (Ndim +dir)*Nproc_+my_rank_;
  int tag2 = (Ndim +dir)*Nproc_+p_recv;

  MPI_Sendrecv( data, size, MPI_DOUBLE, p_send, tag1,
		bin,  size, MPI_DOUBLE, p_recv, tag2,
		MPI_COMM_WORLD, &status );
}
//----------------------------------------------------------
void Communicator::transfer_fw(Field& bin,const Field& data,int size,int dir){
  double abin[size];
  double adata[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_fw(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin.set(i, abin[i]);
}
//----------------------------------------------------------
void Communicator::transfer_bk(Field& bin,const Field& data,int size,int dir){
  double abin[size];
  double adata[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_bk(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin.set(i, abin[i]);
}
//----------------------------------------------------------
void Communicator::transfer_fw(valarray<double>& bin, 
			       const valarray<double>& data, 
			       int size, int dir){
  double abin[size];
  double adata[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_fw(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin[i] = abin[i];
}
//----------------------------------------------------------
void Communicator::transfer_bk(valarray<double>& bin, 
			       const valarray<double>& data,
			       int size,int dir){
  double abin[size];
  double adata[size];
  for(int i=0; i<size; ++i) adata[i] = data[i];
  transfer_bk(abin,adata,size,dir);
  for(int i=0; i<size; ++i) bin[i] = abin[i];
}
//----------------------------------------------------------
void Communicator::send_1to1(double *bin,double *data,int size,
			     int p_to,int p_from,int tag){
  MPI_Status status;
  if(p_to == p_from){
    for(int i=0; i<size; ++i) bin[i] = data[i];
  }else{
   if( my_rank_ == p_from)
     MPI_Send(data, size, MPI_DOUBLE, p_to, tag, MPI_COMM_WORLD );
   if( my_rank_ == p_to)
     MPI_Recv( bin, size, MPI_DOUBLE, p_from, tag, MPI_COMM_WORLD, &status );
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
//----------------------------------------------------------
void Communicator::send_1to1(valarray<double>& bin, 
			     const valarray<double>& data, int size,
			     int p_to, int p_from, int tag){
  if(p_to == p_from){
    bin = data;
  }else{
    double abin[size], adata[size];
    if(my_rank_==p_from) for(int i=0; i<size; ++i) adata[i] = data[i];
    send_1to1(abin,adata,size,p_to,p_from,tag);
    if(my_rank_==p_to)  for(int i=0; i<size; ++i) bin[i] = abin[i];
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
//----------------------------------------------------------
void Communicator::send_1to1(Field& bin, const Field& data, int size,
			     int p_to, int p_from, int tag){
  if(p_to == p_from){
    bin = data;
  }else{
    double abin[size], adata[size];
    if(my_rank_==p_from) for(int i=0; i<size; ++i) adata[i] = data[i];
    send_1to1(abin,adata,size,p_to,p_from,tag);
    if(my_rank_==p_to)  for(int i=0; i<size; ++i) bin.set(i, abin[i]);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
//----------------------------------------------------------
double Communicator::reduce_sum(double& a){
 double a_sum = 0.0;
 MPI_Allreduce(&a, &a_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 return a_sum;
}

//----------------------------------------------------------
void Communicator::sync(){ MPI_Barrier(MPI_COMM_WORLD);}

//----------------------------------------------------------
void Communicator::broadcast(int size, int &data, int sender){
 MPI_Bcast( &data, size, MPI_INT, sender, MPI_COMM_WORLD );
}

//----------------------------------------------------------
void Communicator::broadcast(int size, double &data, int sender){
  MPI_Bcast( &data, size, MPI_DOUBLE, sender, MPI_COMM_WORLD );
}
//----------------------------------------------------------
int Communicator::pprintf(const char* format ...){
  va_list ap;
  int ret = 0;

  va_start(ap,format);
  if(my_rank_==0)
    ret = vfprintf( stdout, format, ap );
  va_end(ap);

  return ret;
}
//**************************************************END*****

