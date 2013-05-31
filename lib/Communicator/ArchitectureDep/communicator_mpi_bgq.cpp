/*!
 * @file communicator_mpi_bgq.cpp
 *
 * @brief Definition of parallel environment Communicator classes, BGQ version
 *
 * Time-stamp: <2013-05-20 23:40:22 cossu>
 *
 */

#include "include/commonPrms.hpp"
#include "Communicator/communicator.hpp"
#include "Communicator/comm_io.hpp"

#include "bgnet.h"
#include <omp.h>

int Communicator::my_rank_;
int Communicator::Nproc_;
int Communicator::ipe_[Ndim_max_];
int Communicator::nd_up_[Ndim_max_];
int Communicator::nd_dn_[Ndim_max_];

static double* pTmpBuffer = NULL;
static int TmpBufferSize = 0;

Communicator* Communicator::instance(){
  static Communicator instance_;
  return &instance_;
}

void Communicator::setup(){
  BGNET_Init();
  my_rank_ = BGNET_Rank();
  Nproc_ = BGNET_Procs();

  int nodeid = my_rank_;

  int NPEx = CommonPrms::instance()->NPEx();
  int NPEy = CommonPrms::instance()->NPEy();
  int NPEz = CommonPrms::instance()->NPEz();
  int NPEt = CommonPrms::instance()->NPEt();

  //Check number of nodes
  if (NPEx*NPEy*NPEz*NPEt != Nproc_) {
    if(my_rank_==0) {
      std::cerr << "Total number of nodes provided in the input file is different from MPI environment ["
		<< Nproc_<<"]\n";
      exit(1);
    }
  }

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

  if (my_rank_ == 0){
    std::cout << "Communicator initialized using MPI with "<< Nproc_ ;
    if (Nproc_==1) 
      std::cout <<" process.\n";
    else
      std::cout <<" processes.\n";
  }

  BGNET_GlobalBarrier();
  CCIO::cout.init(&std::cout);
  CCIO::cerr.init(&std::cerr); 
}

Communicator::~Communicator(){  BGNET_Finalize();}

bool Communicator::primaryNode(){
  if (my_rank_==0) return true;
  return false;
}

//bool Communicator::primaryNode(){ return (my_rank_? true : false);}

int Communicator::nodeid(int x, int y, int z, int t)const{
  CommonPrms* cPar = CommonPrms::instance();
  return x +(cPar->NPEx())*(y +(cPar->NPEy())*(z +(cPar->NPEz())*t));
}

void Communicator::
transfer_fw(double *bin,double *data,int size,int dir) const{
  int p_send = nd_dn_[dir];
  int p_recv = nd_up_[dir];
  
  BGNET_Sendrecv(0,
		 data, size*sizeof(double), p_send,
		 bin, size*sizeof(double), p_recv);
  
}

/*
void Communicator::
transfer_fw_async(double *bin,double *data,int size,int dir) const{
  int p_recv = nd_up_[dir];
  // data is the source pointer 
  // we will have 8 buffers identified by their dir
  int sendID = dir;//in the bw direction shift by 4
  int recvID = sendID + 4; 
  uinst64_t send_offset = 0;
  uinst64_t recv_offset = 0;
  uint64_t size_byte = size*sizeof(double);
  int rcounterID = sendID;

  // prepare the ids for the async comm
  BGNET_SetSendBuffer(data,sendID ,size_byte);
  BGNET_SetRecvBuffer(bin, recvID, size_byte);

  //syncronize after initialization
  BGNET_GlobalBarrier();

  BGNET_Put(sendID,sendID, send_offset, size_byte, p_recv,0,recvID,recv_offset,rcounterID);
}
*/

void Communicator::
transfer_fw(varray_double& bin,const varray_double& data,int dir) const{
  transfer_fw(&bin[0],&(const_cast<varray_double& >(data))[0],
	      bin.size(),dir);
}

void Communicator::
transfer_fw(varray_double& bin,const varray_double& data,
	    const vector_int& index, int dir) const
{
  double* pTmp;
  double* pIn;
  int i,size;
  int p_send = nd_dn_[dir];
  int p_recv = nd_up_[dir];
  int tag1 = dir*Nproc_+my_rank_;
  int tag2 = dir*Nproc_+p_recv;

  size = index.size();

  if(size > TmpBufferSize){
    if(pTmpBuffer != NULL){
      free(pTmpBuffer);
    }
    pTmpBuffer = (double*)malloc(size*sizeof(double));
    TmpBufferSize = size;
  }

  pIn = &(const_cast<varray_double& >(data))[0];
  pTmp = pTmpBuffer;
  for(i=0;i<size;i++){
    pTmp[i] = pIn[index[i]];
  }

  BGNET_Sendrecv(0,pTmp,size*sizeof(double),
		 p_send,&(const_cast<varray_double& >(bin))[0],
		 size*sizeof(double),p_recv);
}

void Communicator::
transfer_bk(double *bin,double *data,int size,int dir) const{
  int p_send = nd_up_[dir];
  int p_recv = nd_dn_[dir];

  BGNET_Sendrecv(0,data,size*sizeof(double),p_send,bin,size*sizeof(double),p_recv);
}

/*
void Communicator::
transfer_bk_async(double *bin,double *data,int size,int dir) const{
  int p_recv = nd_dn_[dir];
  // data is the source pointer 
  // we will have 8 buffers identified by their dir
  int sendID = dir+4;//in the bw direction shift by 4
  int recvID = dir; 
  uinst64_t send_offset = 0;
  uinst64_t recv_offset = 0;
  uint64_t size_byte = size*sizeof(double);
  int rcounterID = sendID;

  // prepare the ids for the async comm
  BGNET_SetSendBuffer(data,sendID ,size_byte);
  BGNET_SetRecvBuffer(bin, recvID, size_byte);

  //syncronize after initialization
  BGNET_GlobalBarrier();

  BGNET_Put(sendID,sendID, send_offset, size_byte, p_recv,0,recvID,recv_offset,rcounterID);
}
*/

void Communicator::
transfer_bk(varray_double& bin,const varray_double& data,int dir) const{
  transfer_bk(&(bin[0]),&(const_cast<varray_double& >(data))[0],
	      bin.size(),dir);
}

void Communicator::
transfer_bk(varray_double& bin,const varray_double& data,
	    const vector_int& index, int dir) const
{
  double* pTmp;
  double* pIn;
  int i,size;
  int p_send = nd_up_[dir];
  int p_recv = nd_dn_[dir];

  size = index.size();

  if(size > TmpBufferSize){
    if(pTmpBuffer != NULL){
      free(pTmpBuffer);
    }
    pTmpBuffer = (double*)malloc(size*sizeof(double));
    TmpBufferSize = size;
  }

  pIn = &(const_cast<varray_double& >(data))[0];
  pTmp = pTmpBuffer;
  for(i=0;i<size;i++){
    pTmp[i] = pIn[index[i]];
  }

  BGNET_Sendrecv(0,pTmp,size*sizeof(double),
		 p_send,&(const_cast<varray_double& >(bin))[0],
		 size*sizeof(double),p_recv);
}

void Communicator::send_1to1(double *bin,double *data,int size,
			     int p_to,int p_from,int tag) const{
  if(p_to == p_from){
    for(int i=0; i<size; ++i) bin[i] = data[i];
  }else{
    if( my_rank_ == p_from)
      BGNET_Send(0,data,size*sizeof(double),p_to);
    if( my_rank_ == p_to)
      BGNET_Recv(0,bin,size*sizeof(double),p_from);
  }
  BGNET_GlobalBarrier();
}

void Communicator::send_1to1(varray_double& bin,
			     const varray_double& data, 
			     int size,int p_to,int p_from,int tag) const{
  if(p_to == p_from)    bin = data;
  else send_1to1(&(bin[0]),&(const_cast<varray_double& >(data))[0],
		 size,p_to,p_from,tag);
  BGNET_GlobalBarrier();
}

double Communicator::reduce_sum(double a) const{
  double a_sum = 0.0;
  BGNET_GlobalSum(&a_sum,&a);
  return a_sum;
}

void Communicator::sync() const{ BGNET_GlobalBarrier();}

void Communicator::broadcast(int size, int &data, int sender) const{
  BGNET_BCast(&data,size,BGNET_COLLECTIVE_INT32,sender,BGNET_COMM_WORLD);
}

void Communicator::broadcast(int size, double &data, int sender) const{
  BGNET_BCast(&data,size,BGNET_COLLECTIVE_DOUBLE,sender,BGNET_COMM_WORLD);
}

int Communicator::reduce_max(double& val,int& idx,int size) const{
  /*! size is the maximum value of idx */
  //  VaId vi = {val, my_rank_*size +idx};
  //  VaId vo;
  double vIn = val;
  int idxIn = my_rank_*size +idx;
  int idxOut;
  BGNET_AllReduceLoc(&vIn,&idxIn,&val,&idxOut,1,
		     BGNET_COLLECTIVE_DOUBLE,BGNET_COLLECTIVE_MAX,
		     BGNET_COMM_WORLD);

  idx = idxOut%size;
  return idxOut/size;
}

int Communicator::reduce_min(double& val,int& idx,int size) const{
  /*! size is the maximum value of idx */
  //  VaId vi = {val, my_rank_*size +idx};
  //  VaId vo;
  double vIn = val;
  int idxIn = my_rank_*size +idx;
  int idxOut;
  BGNET_AllReduceLoc(&vIn,&idxIn,&val,&idxOut,1,
		     BGNET_COLLECTIVE_DOUBLE,BGNET_COLLECTIVE_MIN,
		     BGNET_COMM_WORLD);

  idx = idxOut%size;
  return idxOut/size;
}

