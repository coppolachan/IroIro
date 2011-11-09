//---------------------------------------------------------------------
// sunMat.h
//---------------------------------------------------------------------
#ifndef SUNMAT_INCLUDED
#define SUNMAT_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "include/commonPrms.h"
#endif

#include <iostream>
#include <valarray>

class SUNmat{
private:
  int Nc_;
  std::valarray<double> va_;
public:
  SUNmat(double r=0.0):Nc_(CommonPrms::instance()->Nc()){
    va_.resize(2*Nc_*Nc_,r);
  }
  SUNmat(const SUNmat& m):Nc_(CommonPrms::instance()->Nc()),va_(m.va_){}

  explicit SUNmat(const std::valarray<double>& va)
    :Nc_(CommonPrms::instance()->Nc()),va_(va){}

  const std::valarray<double>& getva() const {return va_;}

  SUNmat& dag();
  SUNmat& unity();
  SUNmat& zero();
  SUNmat& xI();
  SUNmat& reunit();

  SUNmat& operator-();
  
  SUNmat& operator=(const SUNmat&);
  SUNmat& operator=(const double&);

  SUNmat& operator+=(const SUNmat&);
  SUNmat& operator+=(const double&);

  SUNmat& operator-=(const SUNmat&);
  SUNmat& operator-=(const double&);

  SUNmat& operator*=(const SUNmat&);
  SUNmat& operator*=(const double&);

  SUNmat& operator/=(const double&);

  int size() const {return va_.size();}

  double r(int c) const {return va_[2*c  ];}
  double i(int c) const {return va_[2*c+1];}

  double r(int c1,int c2) const {return r(Nc_*c1+c2);}
  double i(int c1,int c2) const {return i(Nc_*c1+c2);}

  void setr(int c, double re){va_[2*c  ] = re;}
  void seti(int c, double im){va_[2*c+1] = im;}

  void setr(int c1,int c2,double re){ setr(Nc_*c1+c2, re);}
  void seti(int c1,int c2,double im){ seti(Nc_*c1+c2, im);}

  void set(int c,double re,double im){
    va_[2*c  ] = re;
    va_[2*c+1] = im;
  }
  void set(int c1,int c2,double re,double im){set(Nc_*c1+c2,re,im);}

  void add(int c, double re, double im){
    va_[2*c  ] += re;
    va_[2*c+1] += im;
  }
  void add(int c1,int c2,double re, double im){ add(Nc_*c1+c2, re, im);}
  
};

inline SUNmat& SUNmat::dag(){
  for(int a = 0; a < Nc_; ++a){
    for(int b = a; b < Nc_; ++b){
      
      int ab = 2*(Nc_*a+b);
      int ba = 2*(Nc_*b+a);

      double re = va_[ab];
      double im = va_[ab+1];

      va_[ab  ] = va_[ba  ];
      va_[ab+1] =-va_[ba+1];

      va_[ba  ] = re;
      va_[ba+1] =-im;
    }
  }
  return *this;
}

inline SUNmat& SUNmat::xI(){
  for(int c = 0; c < va_.size()/2; ++c){
    double tmp = va_[2*c];
    va_[2*c  ] = -va_[2*c+1];
    va_[2*c+1] =  tmp;
  }
  return *this;
}

inline SUNmat& SUNmat::operator-() {
  va_= -va_;
  return *this;
}
inline SUNmat& SUNmat::operator=(const double& rhs){
  va_= rhs;
  return *this;
}
inline SUNmat& SUNmat::operator+=(const SUNmat& rhs){
  va_+= rhs.va_;
  return *this;
}
inline SUNmat& SUNmat::operator+=(const double& rhs){
  va_+= rhs;
  return *this;
}
inline SUNmat& SUNmat::operator-=(const SUNmat& rhs){
  va_-= rhs.va_;
  return *this;
}
inline SUNmat& SUNmat::operator-=(const double& rhs){
  va_-= rhs;
  return *this;
}

inline SUNmat& SUNmat::operator*=(const SUNmat& rhs){
  std::valarray<double> tmp(0.0,2*Nc_*Nc_);

  for(int a = 0; a < Nc_; ++a){
    for(int b = 0; b < Nc_; ++b){
      int ab = 2*(Nc_*a+b);

      for(int c = 0; c < Nc_; ++c){
	int ac = 2*(Nc_*a+c);
	int cb = 2*(Nc_*c+b);
	
	tmp[ab  ]+= va_[ac  ]*rhs.va_[cb  ];
	tmp[ab  ]-= va_[ac+1]*rhs.va_[cb+1];
	tmp[ab+1]+= va_[ac+1]*rhs.va_[cb  ];
	tmp[ab+1]+= va_[ac  ]*rhs.va_[cb+1];
      }
    }
  }
  va_= tmp;
  return *this;
}
inline SUNmat& SUNmat::operator*=(const double& rhs){
  va_*= rhs;
  return *this;
}
inline SUNmat& SUNmat::operator/=(const double& rhs){
  va_ /= rhs;
  return *this;
}

namespace SUNmat_utils{
  SUNmat unity();
  SUNmat zero();

  double ReTr(const SUNmat& m);
  double ImTr(const SUNmat& m);

  const SUNmat dag(const SUNmat& u);
  const SUNmat xI(const SUNmat& u);

  const SUNmat operator+(const SUNmat& m1, const SUNmat& m2);
  const SUNmat operator-(const SUNmat& m1, const SUNmat& m2);
  const SUNmat operator*(const SUNmat& m1, const SUNmat& m2);

  const SUNmat reunit(const SUNmat& m);
  const std::valarray<double> anti_hermite(const SUNmat& m);
}

#endif
