/*!@file realMatrix.hpp
  @brief utility functions for real square matrix
*/
#ifndef REALMATRIX_INCLUDED
#define REALMATRIX_INCLUDED

#include <valarray>

namespace RealMatrix{
  void transpose(std::valarray<double>&,const std::valarray<double>&);
  void invert(std::valarray<double>&,const std::valarray<double>&);
  void trace(double&,const std::valarray<double>&);
}

#endif
