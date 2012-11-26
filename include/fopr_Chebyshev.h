/*! @file fopr_Chebyshev.h
 *  @brief Definition of Chebyshev operator
 */

#ifndef FOPR_CHEBYSHEV_INCLUDED
#define FOPR_CHEBYSHEV_INCLUDED

#include "fopr.h"
#include "Fields/field_expressions.hpp"

class Fopr_Chebyshev :public Fopr_Herm {
private:
  const Fopr_Herm* Op_;
  const int N_;
  mutable Field p_,q_,r_;
public:
  Fopr_Chebyshev(const Fopr_Herm* Op,int N)
    :Op_(Op),N_(N),p_(Op_->fsize()),q_(Op_->fsize()),r_(Op_->fsize()){}
  
  const Field mult(const Field& f) const{
    using namespace FieldExpression;

    q_= f;
    r_= 0.0;

    Field* p = &p_;
    Field* q = &q_;
    Field* r = &r_;
    
    for(int j=0; j<N_; ++j){
      *p = 2.0*Op_->mult(*q) -(*r);
      Field* s = r;
      r = q;
      q = p;
      p = s;
    }
    *q -= Op_->mult(*r);
    return *q;
  }

  const Field mult_dag(const Field& f) const{ return mult(f);}

  double func(double x) const{

    double y = Op_->func(x);

    double p;
    double q = 1.0;
    double r = 0.0;
  
    for(int j=0; j<N_; ++j){  
      p = 2.0*y*q -r;
      r = q;
      q = p;
    }                            
    return q -y*r;
  }

  size_t fsize()const {return Op_->fsize();}
};
#endif
