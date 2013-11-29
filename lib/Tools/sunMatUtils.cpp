//---------------------------------------------------------------------
/*! @file sunMatUtils.cpp
   @brief \f$SU(N)\f$ Matrices linear algebra Utilities
*/ 
//---------------------------------------------------------------------
#include "sunMatUtils.hpp"
#include "Tools/randNum_MP.h"

using namespace std;

namespace SUNmatUtils{

  SUNmat unity(){ return SUNmat().unity();}
  SUNmat zero(){  return SUNmat(); }

  double ReTr(const SUNmat& m){
    double tr = 0.0;
    for(int c=0; c<NC_; ++c) tr += m.r(c,c);
    return tr;
  }
  
  double ImTr(const SUNmat& m){
    double tr = 0.0;
    for(int c=0; c<NC_; ++c) tr += m.i(c,c);
    return tr;
  }

  SUNmat dag(const SUNmat& u){ return SUNmat(u).dag();}
  SUNmat xI(const SUNmat& u){ return SUNmat(u).xI(); }
  SUNmat operator+(const SUNmat& m1,const SUNmat& m2){return SUNmat(m1)+= m2;}
  SUNmat operator-(const SUNmat& m1,const SUNmat& m2){return SUNmat(m1)-= m2;}
  SUNmat operator*(const SUNmat& m1,const SUNmat& m2){return SUNmat(m1)*= m2;}
  SUNmat operator*(const SUNmat& m1,const complex<double>& cp){   
    SUNmat m;
    for(int c=0; c<NC_*NC_; ++c){
      m.setr(c,m1.r(c)*cp.real() -m1.i(c)*cp.imag());
      m.seti(c,m1.r(c)*cp.imag() +m1.i(c)*cp.real());
    }
    return m;
  }

  SUNmat operator*(const SUNmat& m1,double x){ return SUNmat(m1)*= x;  }
  SUNmat operator/(const SUNmat& m1,double x){ return SUNmat(m1)/= x;  }
  SUNmat reunit(const SUNmat& m){ return SUNmat(m).reunit(); }

  SUNmat exponential(const SUNmat& X,int N,int n){
    return N == n ? unity() : exponential(X,N,n+1)*X/n +unity();
  }

  SUNmat anti_hermite_traceless(const SUNmat& m){
    SUNmat out = unity().xI();
    out *= -ImTr(m)/NC_;
    out += anti_hermite(m);
    return out;
  }

  SUNmat anti_hermite(const SUNmat& m){
    SUNmat out(m);
    return out.anti_hermite();
  }

  void SUNprint(const SUNmat& mat) {
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){ 
	std::cout << "("<<mat.r(a,b)<<","<<mat.i(a,b)<<")   ";
      }
      std::cout << "\n";
    }
  }

  /*! @brief Transverse of the outer product of two vectors
    \f[(A^\dagger \circ B)_{ab}^T = A^*_b B_a \f]  
  */
  SUNmat outer_prod_t(const SUNvec& v,const SUNvec& w){
    SUNmat f;
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){
	f.set(a,b,v.r(b)*w.r(a) +v.i(b)*w.i(a),
	          v.r(b)*w.i(a) -v.i(b)*w.r(a));
      }
    }
    return f;
  }

  const SU3mat lambda1(const SU3mat& u){
    SU3mat res(0.0);
    res.set(0, u.r(1,0),u.i(1,0));
    res.set(1, u.r(1,1),u.i(1,1));
    res.set(2, u.r(1,2),u.i(1,2));

    res.set(3, u.r(0,0),u.i(0,0));
    res.set(4, u.r(0,1),u.i(0,1));
    res.set(5, u.r(0,2),u.i(0,2));
    return res;
  }

  const SU3mat lambda2(const SU3mat& u){
    SU3mat res(0.0);
    res.set(0, u.i(1,0),-u.r(1,0));
    res.set(1, u.i(1,1),-u.r(1,1));
    res.set(2, u.i(1,2),-u.r(1,2));

    res.set(3,-u.i(0,0),u.r(0,0));
    res.set(4,-u.i(0,1),u.r(0,1));
    res.set(5,-u.i(0,2),u.r(0,2));
    return res;
  }

  const SU3mat lambda3(const SU3mat& u){
    SU3mat res(0.0);
    res.set(0, u.r(0,0),u.i(0,0));
    res.set(1, u.r(0,1),u.i(0,1));
    res.set(2, u.r(0,2),u.i(0,2));

    res.set(3,-u.r(1,0),-u.i(1,0));
    res.set(4,-u.r(1,1),-u.i(1,1));
    res.set(5,-u.r(1,2),-u.i(1,2));
    return res;
  }

  const SU3mat lambda4(const SU3mat& u){
    SU3mat res(0.0);
    res.set(0, u.r(2,0),u.i(2,0));
    res.set(1, u.r(2,1),u.i(2,1));
    res.set(2, u.r(2,2),u.i(2,2));

    res.set(6, u.r(0,0),u.i(0,0));
    res.set(7, u.r(0,1),u.i(0,1));
    res.set(8, u.r(0,2),u.i(0,2));
    return res;
  }

  const SU3mat lambda5(const SU3mat& u){
    SU3mat res(0.0);
    res.set(0, u.i(2,0),-u.r(2,0));
    res.set(1, u.i(2,1),-u.r(2,1));
    res.set(2, u.i(2,2),-u.r(2,2));

    res.set(6,-u.i(0,0),u.r(0,0));
    res.set(7,-u.i(0,1),u.r(0,1));
    res.set(8,-u.i(0,2),u.r(0,2));
    return res;
  }

  const SU3mat lambda6(const SU3mat& u){
    SU3mat res(0.0);
    res.set(3, u.r(2,0),u.i(2,0));
    res.set(4, u.r(2,1),u.i(2,1));
    res.set(5, u.r(2,2),u.i(2,2));

    res.set(6, u.r(1,0),u.i(1,0));
    res.set(7, u.r(1,1),u.i(1,1));
    res.set(8, u.r(1,2),u.i(1,2));
    return res;
  }

  const SU3mat lambda7(const SU3mat& u){
    SU3mat res(0.0);
    res.set(3, u.i(2,0),-u.r(2,0));
    res.set(4, u.i(2,1),-u.r(2,1));
    res.set(5, u.i(2,2),-u.r(2,2));

    res.set(6, -u.i(1,0),u.r(1,0));
    res.set(7, -u.i(1,1),u.r(1,1));
    res.set(8, -u.i(1,2),u.r(1,2));
    return res;
  }

  const SU3mat lambda8(const SU3mat& u){
    SU3mat res(u);
    res.set(6, -2.0*u.r(2,0),-2.0*u.i(2,0));
    res.set(7, -2.0*u.r(2,1),-2.0*u.i(2,1));
    res.set(8, -2.0*u.r(2,2),-2.0*u.i(2,2));
    res /= sqrt(3.0);
    return res;
  }

  const SU3mat lambda1(){
    SU3mat res(0.0);
    res.setr(1, 1.0);
    res.setr(3, 1.0);
    return res;
  }

  const SU3mat lambda2(){
    SU3mat res(0.0);
    res.seti(1,-1.0);
    res.seti(3, 1.0);
    return res;
  }

  const SU3mat lambda3(){
    SU3mat res(0.0);
    res.setr(0, 1.0);
    res.setr(4,-1.0);
    return res;
  }

  const SU3mat lambda4(){
    SU3mat res(0.0);
    res.setr(2, 1.0);
    res.setr(6, 1.0);
    return res;
  }

  const SU3mat lambda5(){
    SU3mat res(0.0);
    res.seti(2,-1.0); 
    res.seti(6, 1.0);
    return res;
  }

  const SU3mat lambda6(){
    SU3mat res(0.0);
    res.setr(5, 1.0);
    res.setr(7, 1.0);
    return res;
  }

  const SU3mat lambda7(){
    SU3mat res(0.0);
    res.seti(5,-1.0);
    res.seti(7, 1.0);
    return res;
  }

  const SU3mat lambda8(){
    SU3mat res(0.0);
    res.setr(0, 1.0/sqrt(3.0));
    res.setr(4, 1.0/sqrt(3.0));
    res.setr(8,-2.0/sqrt(3.0));
    return res;
  }

  const SU3mat lambdaA(){
    SU3mat res(0.0);
    res.setr(0, 1.0);
    res.setr(8,-1.0);
    return res;
  }
  const SU3mat lambdaB(){
    SU3mat res(0.0);
    res.setr(4, 1.0);
    res.setr(8,-1.0);
    return res;
  }

  #if NC_==3
  // specialization for NC_=3
  template<> 
  const valarray<double> adjoint(const SU3mat& u){

    valarray<double> vt(NADJ_*NADJ_);

    for(int a=0; a<NADJ_; ++a){
      SU3mat lu = lambda_mul[a](u);

      for(int b=0; b<NADJ_; ++b)
	vt[NADJ_*a+b] = ReTr(lu*lambda_mul[b](dag(u)))*0.5;
    }
    return vt;
  }
  /*
  const SU3mat lmd_commutator(int a,int b){
    SU3mat res = lambda_mul[a](lambda_mul[b](unity()));
    res -= lambda_mul[b](lambda_mul[a](unity()));
    return res;
  }
  */
  const SU3mat lmd_commutator(int a,int b){
    return xI(lambda[su3alg[a*NADJ_+b]]())*2.0*su3str[a*NADJ_+b];
  }

  // BLAS style optimization specific functions
  // not working, just containers now

  // Matrix-matrix
  const SUNmat gemm(const SUNmat&, const SUNmat&){}
  // Matrix-(matrix^dag)
  const SUNmat gemmd(const SUNmat&, const SUNmat&){}
  // (Matrix^dag)-matrix
  const SUNmat gemdm(const SUNmat&, const SUNmat&){}
  #else
  // generic code 
  // Matrix-matrix
  const SUNmat gemm(const SUNmat&, const SUNmat&){}
  // Matrix-(matrix^dag)
  const SUNmat gemmd(const SUNmat&, const SUNmat&){}
  // (Matrix^dag)-matrix
  const SUNmat gemdm(const SUNmat&, const SUNmat&){}
  #endif  

}//endof namespace SUNmat_utils
