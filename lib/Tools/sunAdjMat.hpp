/*!@file sunAdjMat.hpp
 * @brief adjoint \f$SU(N)\f$ matrices and linear algebra
 */
#ifndef SUNADJMAT_INCLUDED
#define SUNADJMAT_INCLUDED

#include "include/macros.hpp"
#include "Communicator/comm_io.hpp"
#include "sunMatUtils.hpp"
#include <iostream>
#include <valarray>
#include <assert.h>

template <size_t COLORS = NC_>
class SUNadjMatrix{
private:
  size_t Nadj_;
  std::valarray<double> va_;
public:
  explicit SUNadjMatrix(double r=0.0):Nadj_(COLORS*COLORS-1),va_(r,Nadj_*Nadj_){}

  explicit SUNadjMatrix(const std::valarray<double>& va)
    :Nadj_(COLORS*COLORS-1),va_(va){assert(va.size()==Nadj_*Nadj_);}

  explicit SUNadjMatrix(const SUNmatrix<COLORS>& m)
    :Nadj_(COLORS*COLORS-1),va_(SUNmatUtils::adjoint<COLORS>(m)){}

  SUNadjMatrix(const SUNadjMatrix& m):Nadj_(COLORS*COLORS-1),va_(m.getva()){}

  const std::valarray<double>& getva() const {return va_;}

  SUNadjMatrix& operator-();
  
  SUNadjMatrix& operator=(const SUNadjMatrix&);
  SUNadjMatrix& operator=(const double);

  SUNadjMatrix& operator+=(const SUNadjMatrix&);
  SUNadjMatrix& operator+=(const double);

  SUNadjMatrix& operator-=(const SUNadjMatrix&);
  SUNadjMatrix& operator-=(const double);

  SUNadjMatrix& operator*=(const SUNadjMatrix&);
  SUNadjMatrix& operator*=(const double);

  SUNadjMatrix& operator/=(const double);

  SUNadjMatrix& dag();
  SUNadjMatrix& unity();
  SUNadjMatrix& zero(){ va_= 0.0; }
  SUNadjMatrix& reunit();

  static int size(){return (COLORS*COLORS-1)*(COLORS*COLORS-1);}
  
  double operator[](int a) const {return va_[a];}
  double e(int a1,int a2) const {return va_[Nadj_*a1+a2];}

  void set(int a,double r){va_[a] = r;}
  void set(int a1,int a2,double r){ va_[Nadj_*a1+a2] = r;}

  void add(int a,double r){ va_[a] += r;}
  void add(int a1,int a2,double r){ va_[Nadj_*a1+a2] += r;}
 
  void mult(int a,double r){ va_[a] *= r;}
  void mult(int a1,int a2,double r){ va_[Nadj_*a1+a2] += r;}
};

typedef SUNadjMatrix<3> SU3adjMat;
typedef SUNadjMatrix<> SUNadjMat;

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::unity(){
  va_= 0.0;
  for(int c=0; c<Nadj_; ++c) va_[Nadj_*c+c] = 1.0;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::dag(){
  for(int a=0; a<Nadj_; ++a){
    for(int b=a+1; b<Nadj_; ++b){
      double r = va_[Nadj_*a+b];
      va_[Nadj_*a+b] = va_[Nadj_*b+a];
      va_[Nadj_*b+a] = r;
    }
  }
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator-() {
  va_= -va_;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator=(const double rhs){
  va_= rhs; 
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator=(const SUNadjMatrix& rhs){
  va_= rhs.va_; 
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator+=(const SUNadjMatrix& rhs){
  va_+= rhs.va_;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator+=(const double rhs){
  va_+= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator-=(const SUNadjMatrix& rhs){
  va_-= rhs.va_;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator-=(const double rhs){
  va_-= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator*=(const SUNadjMatrix& rhs){

  std::valarray<double> tmp(0.0,Nadj_*Nadj_);
  for(int a=0; a<Nadj_; ++a){
    for(int b=0; b<Nadj_; ++b){
      int ab = Nadj_*a+b;
      for(int c=0; c<Nadj_; ++c) tmp[ab] += va_[Nadj_*a+c]*rhs.va_[Nadj_*c+b];
    }
  }
  va_= tmp;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator*=(const double rhs){
  va_*= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::operator/=(const double rhs){
  va_/= rhs;
  return *this;
}

// modified Gram-Schmidt method
template <size_t COLORS>
SUNadjMatrix<COLORS>& SUNadjMatrix<COLORS>::reunit(){

  std::valarray<double> u(Nadj_);
  for(int a=0; a<Nadj_; ++a){
    for(int c=0; c<Nadj_; ++c) u[c] = va_[a*Nadj_+c];
    double nrm_i = 1.0/sqrt((u*u).sum());
    for(int c=0; c<Nadj_; ++c) va_[a*Nadj_+c] = u[c]*nrm_i;

    for(int b=a+1; b<Nadj_; ++b){
      double pr = 0.0;
      for(int c=0; c<Nadj_; ++c) pr += va_[a*Nadj_+c]*va_[b*Nadj_+c];
      for(int c=0; c<Nadj_; ++c) va_[b*Nadj_+c] -= pr*va_[a*Nadj_+c];
    }
  }
  return *this;
}

#endif
