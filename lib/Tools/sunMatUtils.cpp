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

  const SUNmat dag(const SUNmat& u){ return SUNmat(u).dag();}
  
  const SUNmat xI(const SUNmat& u){ return SUNmat(u).xI(); }
  
  const SUNmat operator+(const SUNmat& m1,const SUNmat& m2){
    return SUNmat(m1)+= m2;  }

  const SUNmat operator-(const SUNmat& m1,const SUNmat& m2){
    return SUNmat(m1)-= m2;  }

  const SUNmat operator*(const SUNmat& m1,const SUNmat& m2){
    return SUNmat(m1)*= m2;  }

  const SUNmat operator*(const SUNmat& m1,const complex<double>& cp){
    SUNmat m = m1*cp.real();
    return m.xI()*cp.imag();
  }

  const SUNmat operator*(const SUNmat& m1,double x){
    return SUNmat(m1)*= x;  }

  const SUNmat operator/(const SUNmat& m1,double x){
    return SUNmat(m1)/= x;  }

  const SUNmat reunit(const SUNmat& m){ return SUNmat(m).reunit(); }

  /*
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
  */

  /*
  const SUNmat exponential(const SUNmat& X,int N,int n){
    return N == n ? unity() : reunit(exponential(X,N,n+1)*X/n +unity());
  }
  */
  const SUNmat exponential(const SUNmat& X,int N,int n){
    return N == n ? unity() : exponential(X,N,n+1)*X/n +unity();
  }

  const SUNmat anti_hermite_traceless(const SUNmat& m){
    SUNmat out = unity().xI();
    out *= -ImTr(m)/NC_;
    out += anti_hermite(m);
    return out;
  }

  const SUNmat anti_hermite(const SUNmat& m){
    SUNmat out;
    for(int a=0; a<NC_; ++a){
      for(int b=a; b<NC_; ++b){
	double re = 0.5*(m.r(a,b) -m.r(b,a));
	double im = 0.5*(m.i(a,b) +m.i(b,a));
	out.set(a,b, re,im);
	out.set(b,a,-re,im);
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

  /*! @brief Transverse of the outer product of two vectors
    \f[(A^\dagger \circ B)_{ab}^T = A^*_b B_a \f]  
  */
  const SUNmat outer_prod_t(const SUNvec& v,const SUNvec& w){
    SUNmat f;
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){
	f.set(a,b,v.r(b)*w.r(a) +v.i(b)*w.i(a),
	          v.r(b)*w.i(a) -v.i(b)*w.r(a));
      }
    }
    return f;
  }
  /*
  // Creates SUNmat from gausian random number 
  const SUNmat random_mat(const RandNum& rand){
    std::valarray<double> rn(SUNmat::size());
    MPrand::mp_get_gauss(rn,rand);

    return SUNmat(rn).reunit();
  }
  */

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

  template<> 
  const valarray<double> adjoint(const SU3mat& u){
    int dim = 8;
    valarray<double> vt(dim*dim);
    
    for(int a=0; a<dim; ++a){
      SU3mat ua = lambda[a](u);
      for(int b=0; b<dim; ++b){
	SU3mat ub_dag = lambda[b](dag(u));
	ub_dag *= ua;
	vt[dim*a+b] = ReTr(ub_dag)*0.5;
      }
    }
    return vt;
  }

  // BLAS style optimization specific functions
  // not working, just containers now

  #if NC_==3 
  // specialization for NC_=3
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
