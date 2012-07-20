//---------------------------------------------------------------------
/*! @file sunMatUtils.hpp
  @brief \f$SU(N)\f$ Matrices linear algebra Utilities

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef SUNMAT_UTILS_H_
#define SUNMAT_UTILS_H_

#include "sunMat.hpp"
#include "include/common_fields.hpp"

class RandNum;

namespace SUNmatUtils{
  const SUNmat unity();
  const SUNmat zero();
  double ReTr(const SUNmat&);  
  double ImTr(const SUNmat&);
  const SUNmat dag(const SUNmat&);  
  const SUNmat xI(const SUNmat&);  
  const SUNmat operator+(const SUNmat&,const SUNmat&);
  const SUNmat operator-(const SUNmat&,const SUNmat&);
  const SUNmat operator*(const SUNmat&,const SUNmat&);
  const SUNmat operator*(const SUNmat&,double);
  const SUNmat operator/(const SUNmat&,double);
  const SUNmat reunit(const SUNmat&);
  const std::valarray<double> trace_less(const SUNmat&);
  const SUNmat anti_hermite_traceless(const SUNmat&);
  const SUNmat anti_hermite(const SUNmat&);
  const SUNmat outer_prod(const SUNvec& v,const SUNvec& w);
  const SUNmat random_mat(const RandNum& rand);
  const SUNmat exponential(const SUNmat& X,int N,int n=1);

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
