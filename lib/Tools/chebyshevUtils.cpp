#include "chebyshevUtils.hpp"
#include "include/numerical_const.hpp"
#include "Communicator/comm_io.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

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

  /// calculation of T_n(y) for n=0,1,..,N
  double func(double y,vector<double>& Tcb){  
    int N = Tcb.size()-1;
    double p;
    double q = 1.0;
    double r = 0.0;
    Tcb[0] = q;

    for(int j=0; j<N; ++j){  
      p = 2.0*y*q -r;
      r = q;
      q = p;
      Tcb[j+1] = q -y*r;
    }                            
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

  /// getting c'_i's  in f'(x) = sum_{i=0}^{N-1} c'_i*T_i(y)
  /// from c_i's in f(x) = sum_{i=0}^{N} c_i*T_i(y)
  vector<double> dcoeff(const vector<double>& c){
    int N = c.size()-1;
    vector<double> cd(N);

    cd[N-1] = 2.0*double(N)*c[N];
    cd[N-2] = 2.0*double(N-1)*c[N-1];

    for(int j=N-3; j>=0; --j)
      cd[j] = cd[j+2] +2.0*double(j+1)*c[j+1];

    cd[0] *= 0.5;
    return cd;
  }

  double zeros(int l, int n){ return cos((l-0.5)*PI/n);}
  double extrema(int l, int n){  return cos(l*PI/n);  }

  void approx(std::vector<double>& c,int Napprox,ApproxFunc func, 
	      double x0,double x1,void* prms){

    /* Napprox is the order of Chebyshev approximation.
       This function returns the first c.size() coefficients.
       That means the size of c must be set as desired in advance.*/
    
    gsl_cheb_series *cs = gsl_cheb_alloc (Napprox);
    gsl_function F;
    F.function = func;
    F.params = prms;
    gsl_cheb_init(cs,&F,x0,x1);

    c[0] = cs->c[0]*0.5;   //// this is to match our convention
    for(int i=1; i<c.size(); ++i) c[i] = cs->c[i];
    gsl_cheb_free (cs);
  }
  
}
