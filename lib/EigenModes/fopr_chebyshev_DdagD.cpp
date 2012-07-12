/*!
 * @file fopr_chebyshev_DdagD.h
 * @brief Definition of Chebyshev operator
 */
#include "include/fopr_chebyshev_DdagD.h"
#include "Fields/field_expressions.hpp"

const Field Fopr_Chebyshev_DdagD::mult(const Field& f) const{

  using namespace FieldExpression;

  std::vector<Field> dj(3);
  for(int k=0; k< dj.size(); ++k) dj[k].resize(f.size());

  dj[0] = -1.0*f;
  dj[1] = 0.0; 
  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for(int j=Params.Npoly; j>=2; --j){  
    dj[jn] = 2.0*(Params.Fcb1*D_->mult_dag(D_->mult(dj[jp1])) 
		  +Params.Fcb2*dj[jp1]) -dj[jp2];
    jn  = (jn +1)%3;         
    jp1 = (jp1+1)%3;         
    jp2 = (jp2+1)%3;         
  }            
  Field v(f.size());
  v = Params.Fcb1*D_->mult_dag(D_->mult(dj[jp1]))+Params.Fcb2*dj[jp1]-dj[jp2];
  
  return v;
}

double Fopr_Chebyshev_DdagD::func(double x) const{

  std::vector<double> dj(3);
  dj[0] = -1.0;
  dj[1] = 0.0; 
  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for(int j=Params.Npoly; j>=2; --j){  
    dj[jn] = 2.0*(x*x*Params.Fcb1+Params.Fcb2)*dj[jp1] -dj[jp2];

    jn  = (jn +1)%3;         
    jp1 = (jp1+1)%3;         
    jp2 = (jp2+1)%3;         
  }                            
  return (x*x*Params.Fcb1+Params.Fcb2)*dj[jp1] -dj[jp2];              
}
