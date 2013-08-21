//---------------------------------------------------------------------
/*! @file sunVec.hpp
  @brief \f$SU(N)\f$ vectors linear algebra

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef SUNVEC_INCLUDED
#define SUNVEC_INCLUDED

#include "sunMat.hpp"

template <size_t COLORS>
class SUNvector{

private:
  std::valarray<double> va_;
public:
  explicit SUNvector(double r=0.0):va_(r,2*COLORS){}

  explicit SUNvector(const std::valarray<double>& va):va_(va){} 
  
  SUNvector(const SUNvector& v):va_(v.va_){}
  
  const std::valarray<double>& getva() const { return va_;}

  double norm();
  double operator*(const SUNvector&);
  double im_prod(const SUNvector&);
  int nc() const { return COLORS; };
  SUNvector& dag();
  SUNvector& zero();
  SUNvector& xI();

  SUNvector& operator-();

  SUNvector& operator=(double);
  SUNvector& operator=(const std::valarray<double>&);
  SUNvector& operator+=(const SUNvector&);
  SUNvector& operator-=(const SUNvector&);
  SUNvector& operator*=(double);
  SUNvector& operator/=(double);



  int size() const {return 2*COLORS;}

  double r(const int c) const {return va_[2*c  ];}
  double i(const int c) const {return va_[2*c+1];}

  void setr(int c,double re){va_[2*c  ] = re;}
  void seti(int c,double im){va_[2*c+1] = im;}
  void set(int c,double re,double im){
    va_[2*c  ] = re;
    va_[2*c+1] = im;
  }
};

template <size_t COLORS> 
inline double SUNvector<COLORS>::norm(){ 
  return (va_*va_).sum();
}

template <size_t COLORS> 
inline double SUNvector<COLORS>::operator*(const SUNvector& rhs){
  //std::valarray<double> tmp = va_*rhs.va_;
  return (va_*rhs.va_).sum();
}

template <size_t COLORS> 
inline double SUNvector<COLORS>::im_prod(const SUNvector& rhs){
  std::slice re(0,COLORS,2);
  std::slice im(1,COLORS,2);

  std::valarray<double> lhs_im = va_[im];
  std::valarray<double> lhs_re = va_[re];

  return (lhs_re*rhs.va_[im]).sum()-(lhs_im*rhs.va_[re]).sum();
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::dag(){
  for(int c = 0; c < COLORS; ++c) va_[2*c+1] = -va_[2*c+1];
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::zero(){
  va_= 0.0;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::xI(){
  for(int c = 0; c < va_.size()/2; ++c){
    double tmp = va_[2*c];
    va_[2*c  ] = -va_[2*c+1];
    va_[2*c+1] = tmp;
  }
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator-(){
  va_= -va_;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator=(double rhs){
  va_= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator=(const std::valarray<double>& rhs){
  va_.resize(rhs.size());
  va_= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator+=(const SUNvector& rhs){
  va_+= rhs.va_;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator-=(const SUNvector& rhs){
  va_-= rhs.va_;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator*=(double rhs){
  va_*= rhs;
  return *this;
}

template <size_t COLORS> 
inline SUNvector<COLORS>& SUNvector<COLORS>::operator/=(double rhs){
  va_/= rhs;
  return *this;
}

typedef SUNvector<NC_> SUNvec;



namespace SUNvecUtils{
  inline const SUNvec Ix(const SUNvec& u){
    SUNvec tmp;
    for(int c = 0; c < NC_; ++c)
      tmp.set(c, -u.i(c), u.r(c));
    return tmp;
  }

  inline const SUNvec operator+(const SUNvec& v1,const SUNvec& v2){
    return SUNvec(v1)+= v2;
  }
  inline const SUNvec operator-(const SUNvec& v1,const SUNvec& v2){
    return SUNvec(v1)-= v2;
  }
  inline const SUNvec operator*(const SUNvec& v,double r){
    return SUNvec(v)*= r;
  }
  inline const SUNvec operator*(double r, const SUNvec& v){
    return SUNvec(v)*= r;
  }
  inline const SUNvec operator/(const SUNvec& v,double r){
    return SUNvec(v)/= r;
  }

  inline const SUNvec operator*(const SUNmat& m, const SUNvec& v){
    SUNvec tmp;
    
    for(int a = 0; a < NC_; ++a){
      double re = 0.0; 
      double im = 0.0;
      for(int b = 0; b < NC_; ++b){
	re+= m.r(a,b)*v.r(b) -m.i(a,b)*v.i(b);
	im+= m.r(a,b)*v.i(b) +m.i(a,b)*v.r(b);
      }
      tmp.set(a,re,im);
    }
    return tmp;
  }

}

#endif
