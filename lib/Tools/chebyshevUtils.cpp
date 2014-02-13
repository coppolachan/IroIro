#include "signApprox_Chebyshev.hpp"
#include "chebyshevUtils.hpp"

using namespace std;

namespace Chebyshev{
  /// calculation of T_N(y)
  double func(double y,int N){  
    double p;
    double q = 1.0;
    double r = 0.0;
  
    for(int j=0; j<N; ++j){  
      p = 2.0*y*q -r;
      r = q;
      q = p;
    }                            
    return q -y*r;
  }

  /// calculation of f(x) = sum_{i=0}^{N} c_i*T_i(y)
  double series(double y,const vector<double>& c){
    int N = c.size()-1;
    double p;
    double q = c[N];
    double r = 0.0;
  
    for(int j=0; j<N; ++j){  
      p = c[N-j-1] +2.0*y*q -r;
      r = q;
      q = p;
    }                            
    return q -y*r;
  }

  /// getting c'_i in f'(x) = sum_{i=0}^{N-1} c'_i*T_i(y)
  vector<double> dcoeff(const vector<double>& c){
    int N = c.size()-1;
    vector<double> cd(N);

    cd[N-1] = 2.0*double(N)*c[N];
    for(int j=N-1; j>0; --j)
      cd[j-1] = cd[j]+2.0*double(j)*c[j];

    return cd;
  }

}
