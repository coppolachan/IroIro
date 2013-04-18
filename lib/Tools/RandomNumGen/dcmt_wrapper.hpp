/*!
  @file dcmt_wrapper.hpp
  @brief Class to handle the dcmt C library by Matsumoto and Nishimura

 Reference:
 Makoto Matsumoto and Takuji Nishimura,
 "Dynamic Creation of Pseudorandom Number Generators",
 Monte Carlo and Quasi-Monte Carlo Methods 1998,
 Springer, 2000, pp 56--69. 

 Open code here:
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
*/
#ifndef RANDNUM_DCMT_INCLUDED
#define RANDNUM_DCMT_INCLUDED

#include "Tools/randNum.h"
#include "dcmt0.6.1/include/dc.h"

#include<iostream>
#include <string>

class RandNum_DCMT :public RandNum {
private:
  /*
    p should be a Mersenne prime number
    decides the period 2^p -1
    Allowed p values:
    521   607  1279  2203
    2281  3217  4253  4423
    9689  9941 11213 19937
    21701 23209 44497

    w is the output precision 
    31, 32 bit precision only
   */
  enum{w = 32, p = 607};

  mt_struct *mts; /*!< Contains the mt parameters */
  int init(int, int, int);
  double do_rand()const;

  long rand_int31()const;
  
  double rand_double1()const;
  double rand_double2()const;
  double rand_double3()const;
  
  double rand_res53()const;

public:
  // RandNum_DCMT(std::string& file){
  //   loadSeed(file);
  // }
  RandNum_DCMT(int genSeed, int seed, int ID = 0){ init(ID, genSeed, seed);}

  ~RandNum_DCMT();

  void saveSeed(const std::string& file) const;
  void loadSeed(const std::string& file);

};

#endif
