//-------------------------------------------------------------------------
// randNum_MT19937.cpp
//-------------------------------------------------------------------------
#include "randNum_MT19937.h"
#include<iostream>

// add by UEDA @ Nov.2
//-----------------------
#include <string>
#include <fstream>
#include <cstdlib>
#include "Communicator/comm_io.hpp"
//-----------------------

void RandNum_MT19937::init(unsigned long s){
  state_[0]= s & 0xffffffffUL;
  for (int j=1; j<N; ++j) {
    state_[j]= (1812433253UL*(state_[j-1]^(state_[j-1]>>30))+j);
    state_[j]&= 0xffffffffUL;
  }
}

void RandNum_MT19937::init(unsigned long s,
			   unsigned long* key,int key_length){
  init(s);
  int i = 1;
  int j = 0;

  for(int k= N>key_length ? N : key_length; k; --k){
    // non linear
    state_[i] 
      =(state_[i]^((state_[i-1]^(state_[i-1]>>30))*1664525UL))+ key[j] +j; 
    // for WORDSIZE > 32 machines 
    state_[i] &= 0xffffffffUL; 
    
    ++i;
    ++j;
    if(i>=N) { state_[0] = state_[N-1]; i=1; }
    if(j>=key_length) j=0;
  }
  for(int k = N-1; k; --k) {
    // non linear
    state_[i] =(state_[i]^((state_[i-1]^(state_[i-1]>>30))*1566083941UL)) -i; 
    state_[i]&= 0xffffffffUL; // for WORDSIZE > 32 machines 
    i++;
    if(i>=N){state_[0]=state_[N-1]; i=1;}
  }
  // MSB is 1; assuring non-zero initial array 
  state_[0]= 0x80000000UL; 
}

unsigned long RandNum_MT19937::twist(unsigned long u,unsigned long v)const{
  // Period parameters 
  const unsigned long mtrx_a = 0x9908b0dfUL; // constant vector a
  const unsigned long umask = 0x80000000UL;  // most significant w-r bits
  const unsigned long lmask = 0x7fffffffUL;  // least significant r bits 

  unsigned long maxbits = (u & umask)|(v & lmask);
  return (maxbits>>1)^(v&1UL ? mtrx_a : 0UL);
}

void RandNum_MT19937::next_state()const{
 
  unsigned long* p = state_;
  left_= N;
  next_= state_;
  for(int j= N-M+1;--j; ++p) *p = p[M  ]^twist(p[0], p[1]);
  for(int j= M;    --j; ++p) *p = p[M-N]^twist(p[0], p[1]);
  
  *p = p[M-N]^twist(p[0], state_[0]);
}

// generates a random number on [0,0xffffffff]
unsigned long RandNum_MT19937::rand_int32()const{
  if (--left_== 0) next_state();
  unsigned long y = *next_++;

  // Tempering
  y^= (y>>11);
  y^= (y<<7) & 0x9d2c5680UL;
  y^= (y<<15) & 0xefc60000UL;
  y^= (y>>18);
  return y;
}

//generates a random number on [0,0x7fffffff]
long RandNum_MT19937::rand_int31()const{ return (long)(rand_int32()>>1);}

//generates a random number on [0,1] double number
double RandNum_MT19937::rand_double1()const{
  static const double factor = 1.0/4294967295.0;
  return static_cast<double>(rand_int32()*factor);
}

//generates a random number on [0,1) double number
double RandNum_MT19937::rand_double2()const{
  static const double factor = 1.0/4294967296.0;
  return static_cast<double>(rand_int32()*factor);
}

// generates a random number on (0,1) double number
double RandNum_MT19937::rand_double3()const{
  static const double factor = 1.0/4294967296.0;
  return static_cast<double>((rand_int32()+0.5)*factor);
}

// generates a random number on [0,1) with 53-bit resolution
double RandNum_MT19937::rand_res53()const{
  unsigned long a = rand_int32()>>5;
  unsigned long b = rand_int32()>>6;
  static const double factor = 1.0/9007199254740992.0;
  //std::cout<<"randRes53()="<< (a*67108864.0+b)*factor<<std::endl;
  return (a*67108864.0+b)*factor;
}

// add by UEDA @ Nov.2
//--------------------------
// save the seed config
void RandNum_MT19937::saveSeed(std::string& file) {
  CCIO::cout << "Saving Mersenne Twister in file ["
	     << file << "]\n";
  std::ofstream writer(file.c_str());
  for (int i = 0; i < N; ++i) {
    writer << state_[i] << std::endl;
  }
  writer << left_ << std::endl;
}

// load the seed config
void RandNum_MT19937::loadSeed(std::string& file) {
  CCIO::cout << "Loading Mersenne Twister seeds from file ["
	     << file << "]\n";
  
  std::ifstream reader(file.c_str());
  if (reader) {
    unsigned long int n;
    for (int i = 0; i < N; ++i) {
      if (reader >> n) {
        state_[i] = n;
      } else {
        CCIO::cout << "Error: wrong number of seeds in file\n";
        std::abort();
      }
    }
    if (reader >> n) {
      left_ = n;
    } else {
      CCIO::cout << "Error: wrong number of seeds in file\n";
      std::abort();
    }
  } else {
    CCIO::cout << "Error: Seed file ["
	       << file << "] is missing\n";
    std::abort();
  }
  
  next_ = &state_[left_];
}
//--------------------------
