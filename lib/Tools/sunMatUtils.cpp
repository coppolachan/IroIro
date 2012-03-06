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
    return SUNmat(0.0);
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

  const SUNmat anti_hermite_traceless(const SUNmat& m){
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

  const SUNmat anti_hermite(const SUNmat& m){
    std::valarray<double> va(m.getva());
    for(int a=0; a<NC_; ++a){
      for(int b=a; b<NC_; ++b){
	double re = va[2*(NC_*a+b)  ] - va[2*(NC_*b+a)  ];
	double im = va[2*(NC_*a+b)+1] + va[2*(NC_*b+a)+1];
	va[2*(NC_*a+b)  ] =  0.5 * re;
	va[2*(NC_*a+b)+1] =  0.5 * im;
	va[2*(NC_*b+a)  ] = -0.5 * re;
	va[2*(NC_*b+a)+1] =  0.5 * im;
	
      }
    }
    return SUNmat(va);
  }


  void SUNprint(const SUNmat& mat) {
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){ 
	std::cout << "("<<mat.r(a,b)<<","<<mat.i(a,b)<<")   ";
      }
      std::cout << "\n";
    }

  }

/*! @brief Calculates the outer product of two vectors
  \f[(A^\dagger \circ B)_{ab} = A^*_a B_b \f]
 */
  const SUNmat outer_prod(const SUNvec& v,const SUNvec& w){
    
    SUNmat f;
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){
	f.set(a,b,
	      v.r(b)*w.r(a) +v.i(b)*w.i(a),
	      v.r(b)*w.i(a) -v.i(b)*w.r(a));
      }
    }
    return f;
  }


  // BLAS style optimization specific functions
  // not working, just containers now

  #if NC_==3 
  // specialization for NC_=3
  // Matrix-matrix
  const SUNmat gemm(const SUNmat&, const SUNmat&){};
  // Matrix-(matrix^dag)
  const SUNmat gemmd(const SUNmat&, const SUNmat&){};
  // (Matrix^dag)-matrix
  const SUNmat gemdm(const SUNmat&, const SUNmat&){};
  #else
  // generic code 
  // Matrix-matrix
  const SUNmat gemm(const SUNmat&, const SUNmat&){};
  // Matrix-(matrix^dag)
  const SUNmat gemmd(const SUNmat&, const SUNmat&){};
  // (Matrix^dag)-matrix
  const SUNmat gemdm(const SUNmat&, const SUNmat&){};
  #endif  


}//endof namespace SUNmat_utils
