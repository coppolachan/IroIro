//-----------------------------------------------------------------------------
// randNum.h
//----------------------------------------------------------------------------
#ifndef RANDNUM_INCLUDED
#define RANDNUM_INCLUDED

#include <valarray>
#include <math.h>
#include <time.h>
#include <string>



#include "include/numerical_const.hpp"

class RandNum{
private:
  virtual double do_rand() const{return 0;};
public:
  RandNum(){ srand(static_cast<unsigned>(time(NULL)));}
  virtual bool parallel_safe() const{return false;}
  virtual ~RandNum(){}
  
  double get()const{ return do_rand();}


  void get_complex_gauss(double *re, double *im, const double variance = 1.0)const{
    // get the gaussian random number with the variance 1/\sqrt(2)
    // Using the Box-Muller transform
    static double gauss_rand_re = 0;
    static double gauss_rand_im = 0;
    //static double sq2i = 1 / sqrt(2);
    
    //redefine without trigonometric functions...
    double angle = 2*PI*do_rand();
    double rand = 1 -do_rand();
    double factor = sqrt(-log(rand));
    gauss_rand_re = cos(angle)*factor*variance;
    gauss_rand_im = sin(angle)*factor*variance;
    *re = gauss_rand_re;
    *im = gauss_rand_im;
  

  }

  double get_gauss(const double variance = 1.0)const{
    // get the gaussian random number with the variance 1/\sqrt(2)
    static bool has_rand = false;
    static double gauss_rand = 0;
    double temp;
    
    if(has_rand){
      has_rand = false;
      return gauss_rand;
    }else{
      get_complex_gauss(&gauss_rand, &temp, variance);
      /*
      double angle = 2*PI*do_rand();
      double rand = 1 -do_rand();
      double factor = sqrt(-log(rand));
      gauss_rand = cos(angle)*factor;
      */
      has_rand = true;
      //return sin(angle)*factor;
      return temp;
    }
  }
  
  void get(std::valarray<double>& rn) const {
    for(int n=0; n<rn.size(); ++n) rn[n] = do_rand();
  }
  void get_gauss(std::valarray<double>& rn) const{
    for(int n = 0; n < rn.size(); ++n) rn[n] = get_gauss();
  }

  virtual void saveSeed(const std::string&) const = 0;
  virtual void loadSeed(const std::string&) = 0;

};


#endif
