/*! @file fopr_DdagDLinear.h
 *  @brief linear function of DdagD 
 */

#ifndef FOPR_DDAGDLINEAR_INCLUDED
#define FOPR_DDAGDLINEAR_INCLUDED

#include "fopr.h"
#include "Fields/field_expressions.hpp"

class Fopr_DdagDLinear :public Fopr_Herm {
private:
  const Dirac* D_;
  const double a_,b_;
public:
  Fopr_DdagDLinear(const Dirac* D,double a,double b)
    :D_(D),a_(a),b_(b){}
  
  const Field mult(const Field& f) const{
    using namespace FieldExpression;
    
    Field w = D_->mult_dag(D_->mult(f));
    w *= a_;
    w += b_*f;
    return w;
  }
  const Field mult_dag(const Field& f) const{return mult(f);}

  double func(double x)const{return a_*x*x+b_;}
  size_t fsize()const {return D_->fsize();}
};

#endif
