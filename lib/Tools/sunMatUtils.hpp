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

using namespace std;


namespace SUNmatUtils{
  SUNmat unity();
  SUNmat zero();
  double ReTr(const SUNmat&);  
  double ImTr(const SUNmat&);
  const SUNmat dag(const SUNmat&);  
  const SUNmat xI(const SUNmat&);  
  const SUNmat operator+(const SUNmat&, const SUNmat&);
  const SUNmat operator-(const SUNmat&, const SUNmat&);
  const SUNmat operator*(const SUNmat&, const SUNmat&);
  const SUNmat reunit(const SUNmat&);
  const valarray<double> trace_less(const SUNmat&);
  const SUNmat anti_hermite(const SUNmat&);

}//endof namespace SUNmat_utils

#endif
