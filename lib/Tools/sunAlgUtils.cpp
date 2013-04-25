/*!
  @file sunAlgUtils.cpp
  
  @brief Definitions of utility functions which assume the memory alignment
  of the color dof.

  Time-stamp: <2013-04-25 11:48:09 neo>
 */
#include "include/macros.hpp"
#include "sunAlgUtils.hpp"
#include <complex>
#include <cmath>

using namespace std;

namespace SUNmemAling{

  void reunit(double* w){
    complex<double>* pv = (complex<double>*)w;
    complex<double> u[NC_];


    for(int a=0; a<NC_; ++a){
      double nrm = 0.0;
      for(int c=0; c<NC_; ++c){
	u[c] = *(pv+a*NC_+c);
	nrm += abs(u[c])*abs(u[c]);
      }
      double nrm_i = 1.0/sqrt(nrm);
      for(int c=0; c<NC_; ++c) *(pv +a*NC_+c) = u[c]*nrm_i;

      for(int b=a+1; b<NC_; ++b){
	complex<double> z(0.0,0.0);
	for(int c=0; c<NC_; ++c) z += conj(*(pv +a*NC_+c))*(*(pv +b*NC_+c));
	for(int c=0; c<NC_; ++c) *(pv+b*NC_+c) -= z*(*(pv+a*NC_+c));
      }
    }
  }

}
