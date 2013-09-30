//---------------------------------------------------------------------
/*! @file sunMatUtils.hpp
  @brief \f$SU(N)\f$ Matrices linear algebra Utilities
*/ 
//---------------------------------------------------------------------

#ifndef SUNMATUTILS_H_
#define SUNMATUTILS_H_

#include "sunMat.hpp"
#include "sunVec.hpp"
#include "include/errors.hpp"
#include <complex>

class RandNum;

namespace SUNmatUtils{

  // basic functionality 
  SUNmat unity();
  SUNmat zero();
  double ReTr(const SUNmat&);  
  double ImTr(const SUNmat&);
  SUNmat dag(const SUNmat&);  
  SUNmat xI(const SUNmat&);  
  SUNmat operator+(const SUNmat&,const SUNmat&);
  SUNmat operator-(const SUNmat&,const SUNmat&);
  SUNmat operator*(const SUNmat&,const SUNmat&);
  SUNmat operator*(const SUNmat&,double);
  SUNmat operator*(const SUNmat&,const std::complex<double>&);
  SUNmat operator/(const SUNmat&,double);
  SUNmat reunit(const SUNmat&);

  SUNmat anti_hermite_traceless(const SUNmat&);
  SUNmat anti_hermite(const SUNmat&);
  SUNmat outer_prod_t(const SUNvec& v,const SUNvec& w);
  SUNmat exponential(const SUNmat& X,int N,int n=1);

  // multiplication of the Gell-Mann matrices
  const SU3mat lambda1(const SU3mat&);
  const SU3mat lambda2(const SU3mat&);
  const SU3mat lambda3(const SU3mat&);
  const SU3mat lambda4(const SU3mat&);
  const SU3mat lambda5(const SU3mat&);
  const SU3mat lambda6(const SU3mat&);
  const SU3mat lambda7(const SU3mat&);
  const SU3mat lambda8(const SU3mat&);
  
  const SU3mat (*const lambda[])(const SU3mat&) = {
    lambda1,lambda2,lambda3,lambda4,
    lambda5,lambda6,lambda7,lambda8 };

 
 // obtain the adjoint representation
  template<size_t COLORS>
  const std::valarray<double> adjoint(const SUNmatrix<COLORS>& u){
    std::ostringstream msg;
    msg << "implemented only for COLOR=3\n";
    Errors::BaseErr("SU3mat::adjoint", msg);
  }

  template<>
  const std::valarray<double> adjoint(const SUNmatrix<3>& u);

  const SU3mat lmd_commutator(int a,int b);

  //BLAS style multiplications for optimization purposes
  // Matrix-matrix
  const SUNmat gemm(const SUNmat&, const SUNmat&);
  // Matrix-matrix^dag)
  const SUNmat gemmd(const SUNmat&, const SUNmat&);
  // (Matrix^dag)-matrix
  const SUNmat gemdm(const SUNmat&, const SUNmat&);
  //////////////////////////////////////////////////////

  void SUNprint(const SUNmat&);

}//endof namespace SUNmat_utils

#endif
