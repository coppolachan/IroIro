/*!
  @file dcmt_wrapper.cpp
  @brief Class to handle the dcmt C library by Matsumoto and Nishimura

 Reference:
 Makoto Matsumoto and Takuji Nishimura,
 "Dynamic Creation of Pseudorandom Number Generators",
 Monte Carlo and Quasi-Monte Carlo Methods 1998,
 Springer, 2000, pp 56--69. 

 Open code here:
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
*/

#include "dcmt_wrapper.hpp"
#include "include/messages_macros.hpp"
#define BITS_FACTOR32 2.32830643708079737543e-10
#define BITS_FACTOR64 1.11022302462515654042e-16

int RandNum_DCMT::init(int ID, int genSeed, int seed){
  // genSeed 16 bit integer, seed 32 bit
  // Get the independent MT
  _Message(DEBUG_VERB_LEVEL, "Initialization of DC Mersenne Twister\n");
  mts = get_mt_parameter_id_st(w, p, ID, genSeed);
  if (mts == NULL) {
    // print error
  }
  // Initialize RNG mts with 32-bit seed;
  sgenrand_mt(seed, mts);
  return 0;
}

//generates a random number on [0,0x7fffffff]
long RandNum_DCMT::rand_int31()const{ return (long)(genrand_mt(mts)>>1);}

//generates a random number on [0,1] double number
double RandNum_DCMT::rand_double1()const{
  return static_cast<double>(genrand_mt(mts)*BITS_FACTOR32);
}

//generates a random number on [0,1) double number
double RandNum_DCMT::rand_double2()const{
  return static_cast<double>(genrand_mt(mts)*BITS_FACTOR32);
}

// generates a random number on (0,1) double number
double RandNum_DCMT::rand_double3()const{
  return static_cast<double>((genrand_mt(mts)+0.5)*BITS_FACTOR32);
}

// generates a random number on [0,1) with 53-bit resolution
inline double RandNum_DCMT::rand_res53()const{
  return ((genrand_mt(mts)>>5)*67108864.0+(genrand_mt(mts)>>6))*BITS_FACTOR64;
}

double RandNum_DCMT::do_rand()const {
  return rand_res53();
}

RandNum_DCMT::~RandNum_DCMT(){
  free_mt_struct(mts);
}


void RandNum_DCMT::saveSeed(const std::string& file) const{};
void RandNum_DCMT::loadSeed(const std::string& file){};
