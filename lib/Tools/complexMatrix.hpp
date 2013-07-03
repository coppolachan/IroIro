/*!@file complexMatrix.hpp
  @brief utility functions for complex square matrix
*/
#ifndef COMPLEXMATRIX_INCLUDED
#define COMPLEXMATRIX_INCLUDED

#include <valarray>

namespace ComplexMatrix{
  void hermite(std::valarray<double>&,const std::valarray<double>&);
  void transpose(std::valarray<double>&,const std::valarray<double>&);
  void conjugate(std::valarray<double>&,const std::valarray<double>&);
  void invert(std::valarray<double>&,const std::valarray<double>&);
  void trace(double&,double&,const std::valarray<double>&);
}

#endif
