//---------------------------------------------------------------------
/*! @file sunMat.hpp
  @brief \f$SU(N)\f$ Matrices linear algebra

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef SUNMAT_INCLUDED
#define SUNMAT_INCLUDED

#include <iostream>
#include <valarray>
#include <assert.h>

#include "include/macros.hpp"

template <size_t COLORS>
class SUNmatrix{
private:
  std::valarray<double> va_;
public:
  explicit SUNmatrix(double r=0.0):va_(r, 2*COLORS*COLORS){}
  explicit SUNmatrix(const std::valarray<double>& va)
    :va_(va){assert(va.size()==2*COLORS*COLORS);}

  SUNmatrix(const SUNmatrix& m):va_(m.va_){}

  const std::valarray<double>& getva() const {return va_;}

  SUNmatrix& operator-();
  
  SUNmatrix& operator=(const SUNmatrix&);
  SUNmatrix& operator=(const double&);

  SUNmatrix& operator+=(const SUNmatrix&);
  SUNmatrix& operator+=(const double&);

  SUNmatrix& operator-=(const SUNmatrix&);
  SUNmatrix& operator-=(const double&);

  SUNmatrix& operator*=(const SUNmatrix&);
  SUNmatrix& operator*=(const double&);

  SUNmatrix& operator/=(const double&);

  SUNmatrix& dag();
  SUNmatrix& unity();
  SUNmatrix& zero();
  SUNmatrix& xI();
  SUNmatrix& reunit();

  int size() const {return 2*COLORS*COLORS;}

  double r(int c) const {return va_[2*c  ];}
  double i(int c) const {return va_[2*c+1];}

  double r(int c1,int c2) const {return r(COLORS*c1+c2);}
  double i(int c1,int c2) const {return i(COLORS*c1+c2);}

  void setr(int c, double re){va_[2*c  ] = re;}
  void seti(int c, double im){va_[2*c+1] = im;}

  void setr(int c1,int c2,double re){ setr(COLORS*c1+c2, re);}
  void seti(int c1,int c2,double im){ seti(COLORS*c1+c2, im);}

  void set(int c,double re,double im){
    va_[2*c  ] = re;
    va_[2*c+1] = im;
  }
  void set(int c1,int c2,double re,double im){set(COLORS*c1+c2,re,im);}

  void add(int c, double re, double im){
    va_[2*c  ] += re;
    va_[2*c+1] += im;
  }
  void add(int c1,int c2,double re, double im){ add(COLORS*c1+c2, re, im);}
  
};


typedef SUNmatrix<3> SU3mat;
typedef SUNmatrix<NC_> SUNmat;

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::dag(){
  for(int a = 0; a < COLORS; ++a){
    for(int b = a; b < COLORS; ++b){
      
      int ab = 2*(COLORS*a+b);
      int ba = 2*(COLORS*b+a);

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

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::xI(){
  for(int c = 0; c < va_.size()/2; ++c){
    double tmp = va_[2*c];
    va_[2*c  ] = -va_[2*c+1];
    va_[2*c+1] =  tmp;
  }
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator-() {
  va_= -va_;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator=(const double& rhs){
  va_= rhs; 
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator=(const SUNmatrix& rhs){
  va_= rhs.va_; 
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator+=(const SUNmatrix& rhs){
  va_+= rhs.va_;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator+=(const double& rhs){
  va_+= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator-=(const SUNmatrix& rhs){
  va_-= rhs.va_;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator-=(const double& rhs){
  va_-= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator*=(const SUNmatrix& rhs){
  std::valarray<double> tmp(0.0,2*COLORS*COLORS);

  for(int a = 0; a < COLORS; ++a){
    for(int b = 0; b < COLORS; ++b){
      int ab = 2*(COLORS*a+b);

      for(int c = 0; c < COLORS; ++c){
	int ac = 2*(COLORS*a+c);
	int cb = 2*(COLORS*c+b);
	
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

//specialization
//total loop unrolling
template <>
inline SUNmatrix<3>& SUNmatrix<3>::operator*=(const SUNmatrix& rhs){
  std::valarray<double> matrix(18);

  matrix[0]   = va_[0] * rhs.va_[0] - va_[1] * rhs.va_[1] +
                va_[2] * rhs.va_[6] - va_[3] * rhs.va_[7] +
                va_[4] * rhs.va_[12]- va_[5] * rhs.va_[13];
  matrix[1]   = va_[0] * rhs.va_[1] + va_[1] * rhs.va_[0] +
                va_[2] * rhs.va_[7] + va_[3] * rhs.va_[6] +
                va_[4] * rhs.va_[13]+ va_[5] * rhs.va_[12];

  matrix[2]   = va_[0] * rhs.va_[2] - va_[1] * rhs.va_[3] +
                va_[2] * rhs.va_[8] - va_[3] * rhs.va_[9] +
                va_[4] * rhs.va_[14]- va_[5] * rhs.va_[15];
  matrix[3]   = va_[0] * rhs.va_[3] + va_[1] * rhs.va_[2] +
                va_[2] * rhs.va_[9] + va_[3] * rhs.va_[8] +
                va_[4] * rhs.va_[15]+ va_[5] * rhs.va_[14];

  matrix[4]   = va_[0] * rhs.va_[4]  - va_[1] * rhs.va_[5] +
                va_[2] * rhs.va_[10] - va_[3] * rhs.va_[11] +
                va_[4] * rhs.va_[16] - va_[5] * rhs.va_[17];
  matrix[5]   = va_[0] * rhs.va_[5]  + va_[1] * rhs.va_[4] +
                va_[2] * rhs.va_[11] + va_[3] * rhs.va_[10] +
                va_[4] * rhs.va_[17] + va_[5] * rhs.va_[16];



  matrix[6]   = va_[6]  * rhs.va_[0]  - va_[7]  * rhs.va_[1] +
                va_[8]  * rhs.va_[6]  - va_[9]  * rhs.va_[7] +
                va_[10] * rhs.va_[12] - va_[11] * rhs.va_[13];
  matrix[7]   = va_[6]  * rhs.va_[1]  + va_[7]  * rhs.va_[0] +
                va_[8]  * rhs.va_[7]  + va_[9]  * rhs.va_[6] +
                va_[10] * rhs.va_[13] + va_[11] * rhs.va_[12];

  matrix[8]   = va_[6]  * rhs.va_[2]  - va_[7]  * rhs.va_[3] +
                va_[8]  * rhs.va_[8]  - va_[9]  * rhs.va_[9] +
                va_[10] * rhs.va_[14] - va_[11] * rhs.va_[15];
  matrix[9]   = va_[6]  * rhs.va_[3]  + va_[7]  * rhs.va_[2] +
                va_[8]  * rhs.va_[9]  + va_[9]  * rhs.va_[8] +
                va_[10] * rhs.va_[15] + va_[11] * rhs.va_[14];

  matrix[10]  = va_[6]  * rhs.va_[4]  - va_[7]  * rhs.va_[5] +
                va_[8]  * rhs.va_[10] - va_[9]  * rhs.va_[11] +
                va_[10] * rhs.va_[16] - va_[11] * rhs.va_[17];
  matrix[11]  = va_[6]  * rhs.va_[5]  + va_[7]  * rhs.va_[4] +
                va_[8]  * rhs.va_[11] + va_[9]  * rhs.va_[10] +
                va_[10] * rhs.va_[17] + va_[11] * rhs.va_[16];



  matrix[12]  = va_[12] * rhs.va_[0]  - va_[13] * rhs.va_[1] +
                va_[14] * rhs.va_[6]  - va_[15] * rhs.va_[7] +
                va_[16] * rhs.va_[12] - va_[17] * rhs.va_[13];
  matrix[13]  = va_[12] * rhs.va_[1]  + va_[13] * rhs.va_[0] +
                va_[14] * rhs.va_[7]  + va_[15] * rhs.va_[6] +
                va_[16] * rhs.va_[13] + va_[17] * rhs.va_[12];

  matrix[14]  = va_[12] * rhs.va_[2]  - va_[13] * rhs.va_[3] +
                va_[14] * rhs.va_[8]  - va_[15] * rhs.va_[9] +
                va_[16] * rhs.va_[14] - va_[17] * rhs.va_[15];
  matrix[15]  = va_[12] * rhs.va_[3]  + va_[13] * rhs.va_[2] +
                va_[14] * rhs.va_[9]  + va_[15] * rhs.va_[8] +
                va_[16] * rhs.va_[15] + va_[17] * rhs.va_[14];

  matrix[16]  = va_[12] * rhs.va_[4]  - va_[13] * rhs.va_[5] +
                va_[14] * rhs.va_[10] - va_[15] * rhs.va_[11] +
                va_[16] * rhs.va_[16] - va_[17] * rhs.va_[17];
  matrix[17]  = va_[12] * rhs.va_[5]  + va_[13] * rhs.va_[4] +
                va_[14] * rhs.va_[11] + va_[15] * rhs.va_[10] +
                va_[16] * rhs.va_[17] + va_[17] * rhs.va_[16];


  va_= matrix;
  return *this;
}


template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator*=(const double& rhs){
  va_*= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator/=(const double& rhs){
  va_ /= rhs;
  return *this;
}

template <size_t COLORS>
SUNmatrix<COLORS>& SUNmatrix<COLORS>::reunit(){
  double nrm = 0.0;
  for(int c = 0; c < COLORS; ++c) 
    nrm += va_[2*c]*va_[2*c] +va_[2*c+1]*va_[2*c+1];
 
  double nrm_i = 1.0/nrm;
  for(int c = 0; c < COLORS; ++c){
    va_[2*c  ] *= nrm_i;
    va_[2*c+1] *= nrm_i;
  }

  for(int a = 1; a < COLORS; ++a){
    double pr = 0.0;
    double pi = 0.0;

    for(int b = 0; b < a; ++b){
      for(int c = 0; c < COLORS; ++c){
	int ac = a*COLORS+c;
	int bc = b*COLORS+c;
	pr += va_[2*bc]*va_[2*ac  ]+va_[2*bc+1]*va_[2*ac+1];
	pi += va_[2*bc]*va_[2*ac+1]-va_[2*bc+1]*va_[2*ac  ];
      }
      for(int c = 0; c < COLORS; ++c){
	int ac = a*COLORS+c;
	int bc = b*COLORS+c;
	va_[2*ac  ] -= pr*va_[2*bc  ] -pi*va_[2*bc+1];
	va_[2*ac+1] -= pr*va_[2*bc+1] +pi*va_[2*bc  ];
      }
    }
    double nrm = 0.0;
    for(int c = 0; c < COLORS; ++c){
      int ac = a*COLORS+c;
      nrm += va_[2*ac]*va_[2*ac] +va_[2*ac+1]*va_[2*ac+1];
    }
    nrm_i = 1/nrm;
    
    for(int c = 0; c < COLORS; ++c){
      int ac = a*COLORS+c;
      va_[2*ac  ] *= nrm_i;
      va_[2*ac+1] *= nrm_i;
    }
  }
  return *this;
}

#endif
