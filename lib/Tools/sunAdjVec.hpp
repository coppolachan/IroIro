//---------------------------------------------------------------------
/*! @file sunAdjVec.hpp
  @brief \f$SU(N)\f$ad vectors linear algebra
  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef SUNADJVEC_INCLUDED
#define SUNADJVEC_INCLUDED

#include "sunAdjMat.hpp"
#include <valarray>

template <size_t COLORS>
class SUNadjVector{
private:
  std::valarray<double> va_;
public:
  explicit SUNadjVector(double r=0.0):va_(r,COLORS*COLORS-1){}
  explicit SUNadjVector(const std::valarray<double>& va):va_(va){} 
  
  SUNadjVector(const SUNadjVector& v):va_(v.va_){}
  
  const std::valarray<double>& getva() const { return va_;}

  double norm();
  double operator*(const SUNadjVector&);
  SUNadjVector& zero();

  SUNadjVector& operator-();
  SUNadjVector& operator=(double);
  SUNadjVector& operator=(const std::valarray<double>&);
  SUNadjVector& operator+=(const SUNadjVector&);
  SUNadjVector& operator-=(const SUNadjVector&);
  SUNadjVector& operator*=(double);
  SUNadjVector& operator/=(double);

  int size() const {return va_.size();}

  double operator[](const int i) const {return va_[i];}
  void set(int i,double x){va_[i] = x;}
};

template <size_t COLORS> 
inline double SUNadjVector<COLORS>::norm(){ return (va_*va_).sum();}

template <size_t COLORS> 
inline double SUNadjVector<COLORS>::operator*(const SUNadjVector& rhs){
  return (va_*rhs.va_).sum();
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::zero(){
  va_= 0.0;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator-(){
  va_= -va_;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator=(double rhs){
  va_= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator=(const std::valarray<double>& rhs){
  va_.resize(rhs.size());
  va_= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator+=(const SUNadjVector& rhs){
  va_+= rhs.va_;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator-=(const SUNadjVector& rhs){
  va_-= rhs.va_;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator*=(double rhs){
  va_*= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNadjVector<COLORS>& SUNadjVector<COLORS>::operator/=(double rhs){
  va_/= rhs;
  return *this;
}

typedef SUNadjVector<NC_> SUNadjVec;

namespace SUNadjVecUtils{
  inline const SUNadjVec operator+(const SUNadjVec& v1,const SUNadjVec& v2){
    return SUNadjVec(v1)+= v2;
  }
  inline const SUNadjVec operator-(const SUNadjVec& v1,const SUNadjVec& v2){
    return SUNadjVec(v1)-= v2;
  }
  inline const SUNadjVec operator*(const SUNadjVec& v,double r){
    return SUNadjVec(v)*= r;
  }
  inline const SUNadjVec operator*(double r,const SUNadjVec& v){
    return SUNadjVec(v)*= r;
  }
  inline const SUNadjVec operator/(const SUNadjVec& v,double r){
    return SUNadjVec(v)/= r;
  }

  inline const SUNadjVec operator*(const SUNadjMat& m,const SUNadjVec& v){
    SUNadjVec tmp;
    for(int a=0; a<v.size(); ++a){
      double re = 0.0; 
      for(int b=0; b<v.size(); ++b) re+= m.e(a,b)*v[b];
      tmp.set(a,re);
    }
    return tmp;
  }
}

#endif
