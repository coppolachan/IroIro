//-----------------------------------------------------------------------------
/*!
  @file randNum.h
  @brief Declares the abstract class RandNum for random numbers generators
  
*/
//----------------------------------------------------------------------------
#ifndef RANDNUM_INCLUDED
#define RANDNUM_INCLUDED

#include <valarray>
#include <string>



class RandNum{
private:
  // Random number in [0,1)
  virtual double do_rand() const{return 0;}
  
  // Random number in [0,1]
  virtual double do_rand_closed() const{return 0;}

public:
  RandNum(){ srand(static_cast<unsigned>(time(NULL)));}
  virtual bool parallel_safe() const{return false;}
  virtual ~RandNum(){}
  
  double get()const{ return do_rand();}
  double get_gauss(double variance = 1.0)const;
  void get_complex_gauss(double *re,
			 double *im, 
			 double variance = 1.0)const;
  
  void get(std::valarray<double>& rn) const;
  void get_gauss(std::valarray<double>& rn) const;
  
  virtual void saveSeed(const std::string&) const = 0;
  virtual void loadSeed(const std::string&) = 0;

};


#endif
