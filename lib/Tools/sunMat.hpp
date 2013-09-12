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

#define CHECK_TOLERANCE 1e-24

template <size_t COLORS = NC_>
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
  SUNmatrix& operator=(const double);

  SUNmatrix& operator+=(const SUNmatrix&);
  SUNmatrix& operator+=(const double);

  SUNmatrix& operator-=(const SUNmatrix&);
  SUNmatrix& operator-=(const double);

  SUNmatrix& operator*=(const SUNmatrix&);
  SUNmatrix& operator*=(const double);

  SUNmatrix& operator/=(const double);

  SUNmatrix& dag();
  SUNmatrix& star();
  SUNmatrix& unity();
  SUNmatrix& zero() { va_ = 0.0; }
  SUNmatrix& xI();
  SUNmatrix& reunit();
  SUNmatrix& anti_hermite();
  
  bool is_unitary();

  static int size(){return 2*COLORS*COLORS;}
  void print();

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
 
  void mult(int c,double re, double im){
    va_[2*c  ] *= re;
    va_[2*c+1] *= im;
  }
 
  void mult(int c1,int c2,double re, double im){ mult(COLORS*c1+c2, re, im);}
};

typedef SUNmatrix<3> SU3mat;
typedef SUNmatrix<> SUNmat;

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::unity(){
  va_= 0.0;
  for(int c=0; c<COLORS; ++c) va_[2*(COLORS*c+c)] = 1.0;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::dag(){
  for(int a=0; a<COLORS; ++a){
    va_[2*(COLORS*a+a)+1] = -va_[2*(COLORS*a+a)+1];
    for(int b=a+1; b<COLORS; ++b){
      //    for(int b=a; b<COLORS; ++b){
      
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
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::star(){
  for(int a=0; a<COLORS; ++a){
    for(int b=0; b<COLORS; ++b){
      int ab = 2*(COLORS*a+b);
      va_[ab+1] *= -1;
    }
  }
  return *this;
}


template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::anti_hermite(){
  for(int a=0; a<NC_; ++a){
    for(int b=a; b<NC_; ++b){

      int ab = 2*(COLORS*a+b);
      int ba = 2*(COLORS*b+a);

      double re = 0.5*(va_[ab  ]-va_[ba  ]);
      double im = 0.5*(va_[ab+1]+va_[ba+1]);

      va_[ab  ] = re;
      va_[ab+1] = im;

      va_[ba  ] = -re;
      va_[ba+1] = im;
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
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator=(const double rhs){
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
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator+=(const double rhs){
  va_+= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator-=(const SUNmatrix& rhs){
  va_-= rhs.va_;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator-=(const double rhs){
  va_-= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator*=(const SUNmatrix& rhs){

  std::valarray<double> tmp(0.0,2*COLORS*COLORS);
  for(int a=0; a<COLORS; ++a){
    for(int b=0; b<COLORS; ++b){
      int ab = 2*(COLORS*a+b);
      for(int c=0; c<COLORS; ++c){
	int ac = 2*(COLORS*a+c);
	int cb = 2*(COLORS*c+b);
	tmp[ab]  += va_[ac  ]*rhs.va_[cb] -va_[ac+1]*rhs.va_[cb+1];
	tmp[ab+1]+= va_[ac+1]*rhs.va_[cb] +va_[ac  ]*rhs.va_[cb+1];
      }
    }
  }
  va_= tmp;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator*=(double rhs){
  va_*= rhs;
  return *this;
}

template <size_t COLORS>
inline SUNmatrix<COLORS>& SUNmatrix<COLORS>::operator/=(double rhs){
  va_ /= rhs;
  return *this;
}

// modified Gram-Schmidt method
template <size_t COLORS>
SUNmatrix<COLORS>& SUNmatrix<COLORS>::reunit(){
  double nrm_i;
  std::valarray<double> u(2*COLORS);
  for(int a=0; a<COLORS; ++a){
    nrm_i =0.0;
    for(int cc=0; cc<2*COLORS; ++cc)
      {
	//u[cc] = va_[a*2*COLORS +cc];
	nrm_i += va_[a*2*COLORS +cc]*va_[a*2*COLORS +cc];
      }
      nrm_i = 1.0/sqrt(nrm_i);
      for(int cc=0; cc<2*COLORS; cc+=2){
	va_[a*2*COLORS +cc]*=nrm_i;
	va_[a*2*COLORS +cc+1]*=nrm_i;
      }

    for(int b=a+1; b<COLORS; ++b){
      double prr = 0.0;
      double pri = 0.0;
      for(int c=0; c<COLORS; ++c){
	int ac = a*COLORS+c;      
	int bc = b*COLORS+c;
	prr += va_[2*ac]*va_[2*bc  ] +va_[2*ac+1]*va_[2*bc+1];
	pri += va_[2*ac]*va_[2*bc+1] -va_[2*ac+1]*va_[2*bc  ];
      }
      for(int c=0; c<COLORS; ++c){
	int ac = a*COLORS+c;      
	int bc = b*COLORS+c;
	va_[2*bc  ] -= prr*va_[2*ac  ] -pri*va_[2*ac+1];
	va_[2*bc+1] -= prr*va_[2*ac+1] +pri*va_[2*ac  ];
      }
    }
  }
  return *this;
}

template <size_t COLORS>
bool SUNmatrix<COLORS>::is_unitary(){
  double diag_norm_r = 0;
  double diag_norm_i = 0;
  double offdiag_norm_r = 0;
  double offdiag_norm_i = 0;

  std::valarray<double> tmp(0.0,2*COLORS*COLORS);
  for(int a=0; a<COLORS; ++a){
    for(int b=0; b<COLORS; ++b){
      int ab = 2*(COLORS*a+b);
      for(int c=0; c<COLORS; ++c){
	int ac = 2*(COLORS*a+c);
	int bc = 2*(COLORS*b+c);
	tmp[ab]  +=  va_[ac  ]*va_[bc] + va_[ac+1]*va_[bc+1];
	tmp[ab+1]+=  va_[ac+1]*va_[bc] - va_[ac  ]*va_[bc+1];
      }
    }
  }

  for (int a = 0; a < COLORS; a++){
    int aa = 2*(COLORS*a+a);
    diag_norm_r += (tmp[aa]-1.0)*(tmp[aa]-1.0);
    diag_norm_i += tmp[aa+1]*tmp[aa+1];
    
    for (int b = a+1; b < COLORS; b++){
      int ab = 2*(COLORS*a+b);
      int ba = 2*(COLORS*b+a);
      offdiag_norm_r  += tmp[ab]*tmp[ab];
      offdiag_norm_r  += tmp[ba]*tmp[ba];
      offdiag_norm_i  += tmp[ab+1]*tmp[ab+1];
      offdiag_norm_i  += tmp[ba+1]*tmp[ba+1];
    }
  }
 


  if (diag_norm_r > CHECK_TOLERANCE || diag_norm_i > CHECK_TOLERANCE ||
      offdiag_norm_i > CHECK_TOLERANCE || offdiag_norm_i > CHECK_TOLERANCE){
    /*
    CCIO::cout << "diag_norm_r "<< diag_norm_r <<"\n";
    CCIO::cout << "diag_norm_i "<< diag_norm_i <<"\n";
    CCIO::cout << "offdiag_norm_r "<< offdiag_norm_r <<"\n";
    CCIO::cout << "offdiag_norm_i "<< offdiag_norm_i <<"\n";
    */
    return false;

  }
  
  return true;
};


// Utilities 
template <size_t COLORS>
void SUNmatrix<COLORS>::print(){
  for(int a=0; a<COLORS; ++a){
    for(int b=0; b<COLORS; ++b){ 
      std::cout << "("<<this->r(a,b)<<","<<this->i(a,b)<<")   ";
    }
    std::cout << "\n";
  }
  std::cout << "\n"; 
}
#endif
