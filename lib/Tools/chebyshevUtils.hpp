#ifndef CHEBYSHEVUTILS_INCLUDED
#define CHEBYSHEVUTILS_INCLUDED

#include <vector>

namespace Chebyshev{

  typedef double (*ApproxFunc)(double y, void*);
  
  double func(double y,int n);
  double func(double y,std::vector<double>&);
  double series(double y,const std::vector<double>& c);
  double zeros(int l,int n);
  double extrema(int l,int n);
  std::vector<double> dcoeff(const std::vector<double>&);

  void approx(std::vector<double>& c,int Napprox,
	      ApproxFunc,double,double,void* prms=0);
}

#endif
