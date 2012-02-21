//---------------------------------------------------------------------
/*! @file sunMatUtils.cpp
  @brief \f$SU(N)\f$ Matrices linear algebra Utilities

  Class definitions
*/ 
//---------------------------------------------------------------------

#include "sunMatUtils.hpp"

using namespace std;


namespace SUNmatUtils{
  SUNmat unity(){
    SUNmat tmp(0.0);
    for(int c=0; c<NC_; ++c) tmp.setr(c,c,1.0);
    return tmp;
  }

  SUNmat zero(){
    SUNmat tmp(0.0);
    return tmp;
  }

  double ReTr(const SUNmat& m){
    double tr = 0.0;
    for(int c = 0; c < NC_; ++c) tr += m.r(c,c);
    return tr;
  }
  
  double ImTr(const SUNmat& m){
    double tr = 0.0;
    for(int c = 0; c < NC_; ++c) tr += m.i(c,c);
    return tr;
  }

  const SUNmat dag(const SUNmat& u){
    SUNmat tmp(NC_);
    for(int a = 0; a < NC_; a++){
      for(int b = 0; b < NC_; b++){
        tmp.set(a,b, u.r(b,a), -u.i(b,a));
      }
    }
    return tmp;
  }
  
  const SUNmat xI(const SUNmat& u){
    SUNmat tmp(NC_);
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

  const valarray<double> trace_less(const SUNmat& m){
    
    
    double rtr = ReTr(m);
    double itr = ImTr(m);

    rtr /= NC_;
    itr /= NC_;

    valarray<double> va(m.getva());
    for(int c=0; c<NC_; ++c){
      va[2*(NC_*c+c)  ]-= rtr;
      va[2*(NC_*c+c)+1]-= itr;
    }
    return va;
  }

  const SUNmat anti_hermite(const SUNmat& m){
    double trace = ImTr(m);
    trace /= NC_;

    //  valarray<double> va(m.getva());
    SUNmat out;
    for(int a=0; a<NC_; ++a){
      out.add(a,a, 0.0, -trace);
      for(int b=0; b<a; ++b){
	//	int ab = 2*(NC_*a+b);
	//int ba = 2*(NC_*b+a);
	
	out.set(a,b, -out.r(b,a), out.i(b,a));
	out *= 0.5;
	out.set(b,a, -out.r(a,b), out.i(a,b));
	

	/*
	va[ab]-= va[ba];
	va[ab]/= 2.0;
	va[ba] =-va[ab];

	va[ab+1]+= va[ba+1];
	va[ab+1]/= 2.0;
	va[ba+1] = va[ab+1];
	*/
      }
    }
    return out;
  }
}//endof namespace SUNmat_utils
