/*! @file sunAlgUtils.cpp
  @brief Definitions of utility functions which assume the memory alinement
  of the color dof.
 */
#include "include/macros.hpp"
#include "sunAlgUtils.hpp"
#include <complex.h>
#include <math.h>

namespace SUNmemAline{

  void reunit(double* w){
    double _Complex* pv = (double _Complex*)w;
    double _Complex u[NC_];

    for(int a=0; a<NC_; ++a){
      double nrm = 0.0;
      for(int c=0; c<NC_; ++c){
	u[c] = *(pv+a*NC_+c);
	nrm += cabs(u[c])*cabs(u[c]);
      }
      double nrm_i = 1.0/sqrt(nrm);
      for(int c=0; c<NC_; ++c) *(pv +a*NC_+c) = u[c]*nrm_i;

      for(int b=a+1; b<NC_; ++b){
	double _Complex z = 0.0;
	for(int c=0; c<NC_; ++c) z += conj(*(pv +a*NC_+c))*(*(pv +b*NC_+c));
	for(int c=0; c<NC_; ++c) *(pv+b*NC_+c) -= z*(*(pv+a*NC_+c));
      }
    }
  }

}
