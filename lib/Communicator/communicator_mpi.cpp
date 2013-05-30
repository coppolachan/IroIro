/*!
 * @file communicator_mpi.cpp
 *
 * @brief Definition of parallel environment Communicator classes
 *
 * Time-stamp: <2013-05-30 13:00:42 noaki>
 *
 */

#include "communicator.hpp"
#include "commonPrms.hpp"

#include "mpi.h"
#include "comm_io.hpp"

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

  //Check number of nodes
  if (NPEx*NPEy*NPEz*NPEt != Nproc_) {
    if (my_rank_ == 0) {
      std::cerr << "Total number of nodes provided in the input file is different from MPI environment ["
		<< Nproc_<<"]\n";
      exit(1);
      }
  }

  int px = nodeid %NPEx;
  int py = (nodeid/NPEx) %NPEy;
  int pz = (nodeid/(NPEx*NPEy)) %NPEz;
  int pt = nodeid/(NPEx*NPEy*NPEz);

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

  MPI_Barrier(MPI_COMM_WORLD);
  CCIO::cout.init(&std::cout);
  CCIO::cerr.init(&std::cerr);
}

Communicator::~Communicator(){  MPI_Finalize();}

bool Communicator::primaryNode(){
  if (my_rank_==0) return true;
  return false;
}

int Communicator::nodeid(int x, int y, int z, int t)const{
  CommonPrms* cPar = CommonPrms::instance();
  return x +(cPar->NPEx())*(y +(cPar->NPEy())*(z +(cPar->NPEz())*t));
}

void Communicator::allgather(double *bin,double *data,int size)const{
   MPI_Allgather(data,size,MPI_DOUBLE, 
		 bin, size,MPI_DOUBLE,MPI_COMM_WORLD);
}

void Communicator::allgather(varray_double& bin, 
			     const varray_double& data)const{
  allgather(&bin[0],&(const_cast<varray_double& >(data))[0],bin.size());
}

void Communicator::
transfer_fw(double *bin,double *data,int size,int dir)const{
  MPI_Status  status;
  int p_send = nd_dn_[dir];
  int p_recv = nd_up_[dir];

  int tag1 = dir*Nproc_+my_rank_;
  int tag2 = dir*Nproc_+p_recv;

  MPI_Sendrecv(data, size, MPI_DOUBLE, p_send, tag1,
	       bin,  size, MPI_DOUBLE, p_recv, tag2,
	       MPI_COMM_WORLD, &status );
}

void Communicator::
transfer_fw(varray_double& bin,const varray_double& data,int dir)const{
  transfer_fw(&bin[0],&(const_cast<varray_double& >(data))[0],
	      bin.size(),dir);
}

void Communicator::
transfer_fw(varray_double& bin,const varray_double& data,
	    const vector_int& index, int dir)const{
  MPI_Status  status;
  MPI_Datatype subarray;
  MPI_Type_create_indexed_block(index.size(),1,
				&(const_cast<vector_int& >(index))[0],
				MPI_DOUBLE,&subarray); 
  MPI_Type_commit(&subarray);

  int p_send = nd_dn_[dir];
  int p_recv = nd_up_[dir];
  int tag1 = dir*Nproc_+my_rank_;
  int tag2 = dir*Nproc_+p_recv;

  MPI_Sendrecv(&(const_cast<varray_double& >(data))[0],
	       1,subarray,p_send,tag1,
	       &bin[0], index.size(),MPI_DOUBLE,p_recv,tag2,
	       MPI_COMM_WORLD,&status);

  MPI_Type_free (&subarray);
}

void Communicator::
transfer_bk(double *bin,double *data,int size,int dir)const{
  MPI_Status  status;
  int p_send = nd_up_[dir];
  int p_recv = nd_dn_[dir];
  int Ndim = CommonPrms::instance()->Ndim();
  int tag1 = (Ndim +dir)*Nproc_+my_rank_;
  int tag2 = (Ndim +dir)*Nproc_+p_recv;

  MPI_Sendrecv(data, size, MPI_DOUBLE, p_send, tag1,
	       bin,  size, MPI_DOUBLE, p_recv, tag2,
	       MPI_COMM_WORLD, &status );
}

void Communicator::
transfer_bk(varray_double& bin,const varray_double& data,int dir)const{
  transfer_bk(&(bin[0]),&(const_cast<varray_double& >(data))[0],
	      bin.size(),dir);
}

void Communicator::
transfer_bk(varray_double& bin,const varray_double& data,
	    const vector_int& index, int dir)const{
  MPI_Status  status;
  MPI_Datatype subarray;
  int MPIErr;
  MPIErr = MPI_Type_create_indexed_block(index.size(),1,
     				&(const_cast<vector_int& >(index))[0],
     				MPI_DOUBLE,&subarray); 
  MPI_Type_commit(&subarray);

  int p_send = nd_up_[dir];
  int p_recv = nd_dn_[dir];
  int Ndim = CommonPrms::instance()->Ndim();
  int tag1 = (Ndim +dir)*Nproc_+my_rank_;
  int tag2 = (Ndim +dir)*Nproc_+p_recv;

  MPI_Sendrecv(&(const_cast<varray_double& >(data))[0],
	       1,subarray,p_send,tag1,
	       &bin[0], index.size(),MPI_DOUBLE,p_recv,tag2,
	       MPI_COMM_WORLD,&status);

  MPI_Type_free (&subarray);
}

void Communicator::send_1to1(double *bin,double *data,int size,
			     int p_to,int p_from,int tag)const{
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

void Communicator::send_1to1(varray_double& bin,
			     const varray_double& data, 
			     int size,int p_to,int p_from,int tag)const{
  if(p_to == p_from)    bin = data;
  else send_1to1(&(bin[0]),&(const_cast<varray_double& >(data))[0],
		 size,p_to,p_from,tag);
  MPI_Barrier(MPI_COMM_WORLD);
}

double Communicator::reduce_sum(double a)const{
 double a_sum = 0.0;
 MPI_Allreduce(&a, &a_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 return a_sum;
}

void Communicator::sync()const{ MPI_Barrier(MPI_COMM_WORLD);}

void Communicator::broadcast(int size, int &data, int sender)const{
 MPI_Bcast( &data, size, MPI_INT, sender, MPI_COMM_WORLD );
}

void Communicator::broadcast(int size, double &data, int sender)const{
  MPI_Bcast( &data, size, MPI_DOUBLE, sender, MPI_COMM_WORLD );
}

int Communicator::reduce_max(double& val,int& idx,int size)const{
  /*! size is the maximum value of idx */
  VaId vi = {val, my_rank_*size +idx};
  VaId vo;
  MPI_Reduce(&vi,&vo,1,MPI_DOUBLE_INT,MPI_MAXLOC,0,MPI_COMM_WORLD);
  val = vo.value;
  idx = vo.index%size;
  return vo.index/size;
}

int Communicator::reduce_min(double& val,int& idx,int size)const{
  /*! size is the maximum value of idx */
  
  VaId vi = {val, my_rank_*size +idx};
  VaId vo;
  MPI_Reduce(&vi,&vo,1,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
  val = vo.value;
  idx = vo.index%size;
  return vo.index/size;
}




