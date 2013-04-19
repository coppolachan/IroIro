/*!
  @file randNum.cpp

 */

#include <math.h>

#include "randNum.h"
#include "include/numerical_const.hpp"

void RandNum::get_complex_gauss(double *re, double *im, const double variance)const{
  // get the gaussian random number with the variance 1/\sqrt(2)
  // Using the Box-Muller transform
  static double gauss_rand_re = 0;
  static double gauss_rand_im = 0;
  
  //redefine without trigonometric functions...
  /*
    double angle = 2*PI*do_rand();
    double rand = 1 -do_rand();
    double factor = sqrt(-log(rand));
    gauss_rand_re = cos(angle)*factor*variance;
    gauss_rand_im = sin(angle)*factor*variance;
    *re = gauss_rand_re;
    *im = gauss_rand_im;
    */
  
  // Polar form
  double x1, x2, w;
  do {
    x1 = 2.0 * do_rand_closed()- 1.0;
    x2 = 2.0 * do_rand_closed()- 1.0;
      w = x1 * x1 + x2 * x2;
  }while ( w >= 1.0 );
  
  w = sqrt( (-2.0 * log( w ) ) / w );
  *re = x1 * w * variance;
  *im = x2 * w * variance;
}

double RandNum::get_gauss(const double variance)const{
  // get the gaussian random number with the variance 1/\sqrt(2)
  static bool has_rand = false;
  static double gauss_rand = 0;
  double temp;
  
  if(has_rand){
    has_rand = false;
    return gauss_rand;
  }else{
    get_complex_gauss(&gauss_rand, &temp, variance);
    has_rand = true;
    return temp;
  }
}

void RandNum::get(std::valarray<double>& rn) const  {
  for(int n=0; n<rn.size(); ++n) rn[n] = do_rand();
}
void RandNum::get_gauss(std::valarray<double>& rn) const{
  if (rn.size() & 1)
    for(int n = 0; n < rn.size(); ++n)  rn[n] = get_gauss();
  else
    for(int n = 0; n < rn.size(); n+=2){
      get_complex_gauss(&rn[n], &rn[n+1]);
    }
}


