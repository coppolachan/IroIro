//-----------------------------------------------------------------------------
// communicator.h
//-----------------------------------------------------------------------------
#ifndef COMMUNICATOR_INCLUDED
#define COMMUNICATOR_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "include/commonPrms.h"
#endif

#include <valarray>
#include <stdio.h>
#include <cstdarg>

class Field;

class Communicator{
private:
  static int my_rank_;        // rank of current process
  static int Nproc_;          // number of processes

  enum{Ndim_max_ = 4}; // enum hack

  static int ipe_[Ndim_max_];
  static int nd_up_[Ndim_max_];
  static int nd_dn_[Ndim_max_];

  static void setup();
  Communicator(){
    setup(); // initializers and setup are called only once
  }
  Communicator(const Communicator&){}

public:
  ~Communicator();
  static Communicator* instance();

  static int size(){return Nproc_;}
  static int id(){return my_rank_;}
  static int nodeid(){return my_rank_;}
  static int nodeid(int x, int y, int z, int t);
  static int ipe(int dir){ return ipe_[dir];}
  static double reduce_sum(double&);
  static void sync();
  static void broadcast(int size, int &data, int sender);
  static void broadcast(int size, double &data, int sender);

  static void transfer_fw(double *bin, double *data,int size,int dir);
  static void transfer_bk(double *bin, double *data,int size,int dir);

  static void transfer_fw(std::valarray<double>& bin, 
			  const std::valarray<double>& data,int size,int dir);
  static void transfer_bk(std::valarray<double>& bin, 
			  const std::valarray<double>& data,int size,int dir);

  static void transfer_fw(Field& bin, const Field& data,int size,int dir);
  static void transfer_bk(Field& bin, const Field& data,int size,int dir);

  static void send_1to1(double *bin, double *data,int size,
			int p_to, int p_from, int tag);

  static void send_1to1(std::valarray<double>& bin, 
			const std::valarray<double>& data,int size,
			int p_to, int p_from, int tag);

  static void send_1to1(Field& bin, const Field& data,int size,
			int p_to, int p_from, int tag);

  static double get_time(){return 0.0;};

  static int pprintf(const char* ...);

  static bool primaryNode();

};

namespace CommunicatorItems{

 inline int pprintf(const char* format ...){

  va_list ap;
  int ret = 0;

  va_start(ap,format);
  if(Communicator::nodeid()==0) ret = vfprintf(stdout, format, ap);
  va_end(ap);

  return ret;
 }

}
#endif
