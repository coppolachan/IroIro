//---------------------------------------------------------------------
// sunVec.h
//---------------------------------------------------------------------
#ifndef SUNVEC_INCLUDED
#define SUNVEC_INCLUDED

#ifndef SUNMAT_INCLUDED
#include "sunMat.h"
#endif

#include <valarray>

class SUNvec{

private:
  int Nc_;  
  std::valarray<double> va_;
public:
  SUNvec(double r=0.0):Nc_(CommonPrms::instance()->Nc()){
    va_.resize(2*Nc_,r);
  }
  SUNvec(const SUNvec& v):Nc_(CommonPrms::instance()->Nc()){
    va_.resize(v.size());
    va_= v.va_;
  }
  explicit SUNvec(const std::valarray<double>& va)
    :Nc_(CommonPrms::instance()->Nc()){
    va_.resize(va.size());
    va_= va; 
  }
  const std::valarray<double>& getva() const { return va_;}

  double norm();
  double operator*(const SUNvec&);
  int nc() const { return Nc_; };
  SUNvec& dag();
  SUNvec& zero();
  SUNvec& xI();

  SUNvec& operator-();

  SUNvec& operator=(const double&);
  SUNvec& operator=(const std::valarray<double>&);
  //  SUNvec& operator=(const std::complex<double>&);

  SUNvec& operator+=(const SUNvec&);
  SUNvec& operator-=(const SUNvec&);
  SUNvec& operator*=(const double&);
  SUNvec& operator/=(const double&);
  //  SUNvec& operator*=(const std::complex<double>&);
  //  SUNvec& operator/=(const std::complex<double>&);

  /*
  SUNvec& operator=(const MVmult& mv){
    sun_mat_vec(this, mv.m, mv.v);
    return *this;
  }
  */

  int size() const {return va_.size();}

  double r(const int c) const {return va_[2*c  ];}
  double i(const int c) const {return va_[2*c+1];}

  void setr(const int c, const double re){va_[2*c  ] = re;}
  void seti(const int c, const double im){va_[2*c+1] = im;}
  void set(const int c, const double re, const double im){
    va_[2*c  ] = re;
    va_[2*c+1] = im;
  }
};

inline double SUNvec::norm(){ 
  std::valarray<double> tmp = va_*va_;
  return tmp.sum();
}
inline double SUNvec::operator*(const SUNvec& rhs){
  std::valarray<double> tmp = va_*rhs.va_;
  return tmp.sum();
}
inline SUNvec& SUNvec::dag(){
  for(int c = 0; c < Nc_; ++c) va_[2*c+1] = -va_[2*c+1];
  return *this;
}
inline SUNvec& SUNvec::zero(){
  va_= 0.0;
  return *this;
}
inline SUNvec& SUNvec::xI(){
  for(int c = 0; c < va_.size()/2; ++c){
    double tmp = va_[2*c];
    va_[2*c  ] = -va_[2*c+1];
    va_[2*c+1] = tmp;
  }
  return *this;
}
inline SUNvec& SUNvec::operator-(){
  va_= -va_;
  return *this;
}
inline SUNvec& SUNvec::operator=(const double& rhs){
  va_= rhs;
  return *this;
}
inline SUNvec& SUNvec::operator=(const std::valarray<double>& rhs){
  va_.resize(rhs.size());
  va_= rhs;
  return *this;
}
inline SUNvec& SUNvec::operator+=(const SUNvec& rhs){
  va_+= rhs.va_;
  return *this;
}
inline SUNvec& SUNvec::operator-=(const SUNvec& rhs){
  va_-= rhs.va_;
  return *this;
}
inline SUNvec& SUNvec::operator*=(const double& rhs){
  va_*= rhs;
  return *this;
}
/*
inline SUNvec& SUNvec::operator*=(const std::complex<double>& rhs){
  std::valarray<double> tmp = va_;
  for(int c = 0; c < va_.size()/2; ++c){
    va_[2*c  ] = (tmp[2*c]*rhs.real() -tmp[2*c+1]*rhs.imag());
    va_[2*c+1] = (tmp[2*c]*rhs.imag() +tmp[2*c+1]*rhs.real());
  }
  return *this;
}
*/
inline SUNvec& SUNvec::operator/=(const double& rhs){
  va_/= rhs;
  return *this;
}
/*
inline SUNvec& SUNvec::operator/=(const std::complex<double>& rhs){
  std::valarray<double> tmp = va_;
  for(int c = 0; c < va_.size()/2; ++c){
    va_[2*c  ] = ( tmp[2*c]*rhs.real() +tmp[2*c+1]*rhs.imag())/abs(rhs);
    va_[2*c+1] = (-tmp[2*c]*rhs.imag() -tmp[2*c+1]*rhs.real())/abs(rhs);
  }
  return *this;
}
*/

//using namespace SUNmat_utils;

namespace SUNvec_utils{
  inline const SUNvec Ix(const SUNvec& u){
    int Nc_ = u.nc();
    SUNvec tmp(Nc_);
    for(int c = 0; c < u.size()/2; ++c)
      tmp.set(c, -u.i(c), u.r(c));
    return tmp;
  }

  inline const SUNvec operator+(const SUNvec& v1, const SUNvec& v2){
    return SUNvec(v1)+= v2;
  }
  inline const SUNvec operator-(const SUNvec& v1, const SUNvec& v2){
    return SUNvec(v1)-= v2;
  }
  inline const SUNvec operator*(const SUNvec& v, const double& r){
    return SUNvec(v)*= r;
  }
  inline const SUNvec operator*(const double& r, const SUNvec& v){
    return SUNvec(v)*= r;
  }
  inline const SUNvec operator/(const SUNvec& v, const double& r){
    return SUNvec(v)/= r;
  }

  inline const SUNvec operator*(const SUNmat& m, const SUNvec& v){
    int Nc_ = v.size()/2;
    SUNvec tmp(v.size()); 
    
    for(int a = 0; a < Nc_; ++a){
      double re = 0.0; 
      double im = 0.0;
      for(int b = 0; b < Nc_; ++b){
	re+= m.r(a,b)*v.r(b) -m.i(a,b)*v.i(b);
	im+= m.r(a,b)*v.i(b) +m.i(a,b)*v.r(b);
      }
      tmp.set(a,re,im);
    }
    return tmp;
  }

}

#endif
