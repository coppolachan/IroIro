/*! @file sunAdjMatUtils.cpp
 *  @brief adjoint \f$SU(N)\f$ Matrices linear algebra Utilities
*/ 
#include "sunAdjMatUtils.hpp"

using namespace std;

namespace SUNadjMatUtils{

  const SUNadjMat unity(){ return SUNadjMat().unity();}
  const SUNadjMat zero(){  return SUNadjMat(); }

  double Tr(const SUNadjMat& m){
    int dim = sqrt(double(m.size()));
    double tr = 0.0;
    for(int a=0; a<dim; ++a) tr += m.e(a,a);
    return tr;
  }
}
