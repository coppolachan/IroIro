//-------------------------------------------------------------------------
// randNum_MT19937.h
//-------------------------------------------------------------------------
#ifndef RANDNUM_MT19937_INCLUDED
#define RANDNUM_MT19937_INCLUDED

#include<iostream>
#include <string>

#include "Tools/randNum.h"
#include "Communicator/comm_io.hpp"
class RandNum_MT19937 :public RandNum {
private:
  enum{N=624, M=397};

  mutable int left_;
  mutable unsigned long state_[N];
  mutable unsigned long* next_;

  unsigned long rand_int32()const;
  unsigned long twist(unsigned long u, unsigned long v)const;

  void init(unsigned long s);
  void init(unsigned long s, unsigned long* key, int key_length);
  void next_state()const;
  long rand_int31()const;

  double rand_double1()const;
  double rand_double2()const;
  double rand_double3()const;

  double rand_res53()const;
  double do_rand()const{ return rand_res53();}
  double do_rand_closed()const{ return rand_double1();}
public:
  RandNum_MT19937(unsigned long s = 5489UL):left_(1){init(s);}
  RandNum_MT19937(unsigned long* key,int key_length)
    :left_(1){init(19650218UL,key,key_length); }
  RandNum_MT19937(std::string& file){loadSeed(file);}

  bool parallel_safe() const{return false;}

  ~RandNum_MT19937(){};

  void saveSeed(const std::string& file) const;
  void loadSeed(const std::string& file);

};

#endif
