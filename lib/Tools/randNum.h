//-----------------------------------------------------------------------------
// randNum.h
//----------------------------------------------------------------------------
#ifndef RANDNUM_INCLUDED
#define RANDNUM_INCLUDED

#include <valarray>
#include <math.h>
#include <time.h>


static const double PI = 6.0*asin(0.5);
//static const double PI = 3.141592653589793;


class RandNum{
private:
  virtual double do_rand() const{}
public:
  RandNum(){ srand(static_cast<unsigned>(time(NULL)));}
  virtual ~RandNum(){}
  
  double get()const{ return do_rand();}
  double get_gauss()const{
    // get the gaussian random number with the variance 1/\sqrt(2)
    static bool has_rand = false;
    static double gauss_rand = 0;
    //static double sq2i = 1 / sqrt(2);
    
    if(has_rand){
      has_rand = false;
      return gauss_rand;
    }else{
      double angle = 2*PI*do_rand();
      double rand = 1 -do_rand();
      double factor = sqrt(-log(rand));
      gauss_rand = cos(angle)*factor;
      has_rand = true;
      return sin(angle)*factor;
    }
  }
  
  void get(std::valarray<double>& rn) const {
    for(int n=0; n<rn.size(); ++n) rn[n] = do_rand();
  }
  void get_gauss(std::valarray<double>& rn) const{
    for(int n = 0; n < rn.size(); ++n) rn[n] = get_gauss();
  }
};


#endif
