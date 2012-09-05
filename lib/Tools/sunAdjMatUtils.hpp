/*! @file sunAdjMatUtils.hpp
 *  @brief adjoint \f$SU(N)\f$ Matrices linear algebra Utilities
*/ 
#ifndef SUNADJMATUTILS_H_
#define SUNADJMATUTILS_H_

#include "sunAdjMat.hpp"

namespace SUNadjMatUtils{

  const SUNadjMat unity();
  const SUNadjMat zero();

  double Tr(const SUNadjMat&);
}
#endif
