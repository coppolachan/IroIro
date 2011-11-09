//----------------------------------------------------------------------
//   sunMat.cpp
//----------------------------------------------------------------------
#include "sunMat.h"
using namespace std;

SUNmat& SUNmat::reunit(){
  double nrm = 0.0;
  for(int c = 0; c < Nc_; ++c) 
    nrm += va_[2*c]*va_[2*c] +va_[2*c+1]*va_[2*c+1];
 
  double nrm_i = 1.0/nrm;
  for(int c = 0; c < Nc_; ++c){
    va_[2*c  ] *= nrm_i;
    va_[2*c+1] *= nrm_i;
  }

  for(int a = 1; a < Nc_; ++a){
    double pr = 0.0;
    double pi = 0.0;

    for(int b = 0; b < a; ++b){
      for(int c = 0; c < Nc_; ++c){
	int ac = a*Nc_+c;
	int bc = b*Nc_+c;
	pr += va_[2*bc]*va_[2*ac  ]+va_[2*bc+1]*va_[2*ac+1];
	pi += va_[2*bc]*va_[2*ac+1]-va_[2*bc+1]*va_[2*ac  ];
      }
      for(int c = 0; c < Nc_; ++c){
	int ac = a*Nc_+c;
	int bc = b*Nc_+c;
	va_[2*ac  ] -= pr*va_[2*bc  ] -pi*va_[2*bc+1];
	va_[2*ac+1] -= pr*va_[2*bc+1] +pi*va_[2*bc  ];
      }
    }
    double nrm = 0.0;
    for(int c = 0; c < Nc_; ++c){
      int ac = a*Nc_+c;
      nrm += va_[2*ac]*va_[2*ac] +va_[2*ac+1]*va_[2*ac+1];
    }
    nrm_i = 1/nrm;
    
    for(int c = 0; c < Nc_; ++c){
      int ac = a*Nc_+c;
      va_[2*ac  ] *= nrm_i;
      va_[2*ac+1] *= nrm_i;
    }
  }
  return *this;
}

namespace SUNmat_utils{
  SUNmat unity(){
    int Nc = CommonPrms::instance()->Nc();
    SUNmat tmp(0.0);
    for(int c=0; c<Nc; ++c) tmp.setr(c,c,1.0);
    return tmp;
  }

  SUNmat zero(){
    SUNmat tmp(0.0);
    return tmp;
  }

  double ReTr(const SUNmat& m){
    int Nc = CommonPrms::instance()->Nc();
    double tr = 0.0;
    for(int c = 0; c < Nc; ++c) tr += m.r(c,c);
    return tr;
  }
  
  double ImTr(const SUNmat& m){
    int Nc = CommonPrms::instance()->Nc();
    double tr = 0.0;
    for(int c = 0; c < Nc; ++c) tr += m.i(c,c);
    return tr;
  }

  const SUNmat dag(const SUNmat& u){
    int Nc = CommonPrms::instance()->Nc();
    SUNmat tmp(Nc);
    for(int a = 0; a < Nc; a++){
      for(int b = 0; b < Nc; b++){
        tmp.set(a,b, u.r(b,a), -u.i(b,a));
      }
    }
    return tmp;
  }
  
  const SUNmat xI(const SUNmat& u){
    int Nc = CommonPrms::instance()->Nc();
    SUNmat tmp(Nc);
    for(int c = 0; c < u.size()/2; ++c)
      tmp.set(c, -u.i(c), u.r(c));
    return tmp;
  }
  
  const SUNmat operator+(const SUNmat& m1, const SUNmat& m2){
    return SUNmat(m1)+= m2;
  }

  const SUNmat operator-(const SUNmat& m1, const SUNmat& m2){
    return SUNmat(m1)-= m2;
  }

  const SUNmat operator*(const SUNmat& m1, const SUNmat& m2){
    return SUNmat(m1)*= m2;
  }

  const SUNmat reunit(const SUNmat& m){
    SUNmat tmp = m;
    return tmp.reunit();
  }

  const valarray<double> anti_hermite(const SUNmat& m){

    int Nc = CommonPrms::instance()->Nc();
    double trace = ImTr(m);
    trace /= Nc;

    valarray<double> va(m.getva());
    for (int a=0; a<Nc; ++a) {
      va[2*(Nc*a+a)  ]= 0.0;
      va[2*(Nc*a+a)+1]-= trace;
      
      for (int b=0; b<a; ++b) {
	int ab = 2*(Nc*a+b);
	int ba = 2*(Nc*b+a);

	va[ab]-= va[ba];
	va[ab]/= 2.0;
	va[ba] =-va[ab];

	va[ab+1]+= va[ba+1];
	va[ab+1]/= 2.0;
	va[ba+1] = va[ab+1];
      }
    }
    return va;
  }
}

