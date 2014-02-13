#ifndef CHEBYSHEVUTILS_INCLUDED
#define CHEBYSHEVUTILS_INCLUDED

#include <vector>

namespace Chebyshev{

  double func(double y,int N);
  double series(double y,std::vector<double>& c);
  std::vector<double> dcoeff(const std::vector<double>&);
  
}

#endif
