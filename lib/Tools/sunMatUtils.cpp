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

    SUNmat out = m;

    for(int a=0; a<NC_; ++a){
      out.setr(a,a, 0.0);
      out.add(a,a, 0.0, -trace);
       for(int b=0; b<a; ++b){
	out.add(a,b, -out.r(b,a), out.i(b,a));
	out.mult(a,b,0.5,0.5);
	out.set(b,a, -out.r(a,b), out.i(a,b));
      }
      
    }
    return out;
  }

  void SUNprint(const SUNmat& mat) {
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){ 
	std::cout << "("<<mat.r(a,b)<<","<<mat.i(a,b)<<")   ";
      }
      std::cout << "\n";
    }

  }




}//endof namespace SUNmat_utils
