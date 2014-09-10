/*! @file fopr_Chebyshev.h
 *  @brief Definition of Chebyshev operator
 */

#ifndef FOPR_CHEBYSHEV_INCLUDED
#define FOPR_CHEBYSHEV_INCLUDED

#include "fopr.h"
#include "Fields/field_expressions.hpp"
#include "pugi_interface.h"
#include "Tools/chebyshevUtils.hpp"

class Fopr_Chebyshev :public Fopr_Herm {
private:
  const Fopr_Herm* Op_;
  int N_;
  mutable Field p_,q_,r_;
  std::vector<double> c_;
public:
  /// to obtain T_N(x)
  Fopr_Chebyshev(int N,const Fopr_Herm* Op)
    :N_(N),Op_(Op),p_(Op_->fsize()),q_(Op_->fsize()),r_(Op_->fsize()),
     c_(N+1,0.0){
    c_[N] = 1.0; 
  }
  /// to obtain sum_i c_i T_i(x)
  Fopr_Chebyshev(const Fopr_Herm* Op,const std::vector<double>& c)
    :N_(c.size()-1),Op_(Op),
     p_(Op_->fsize()),q_(Op_->fsize()),r_(Op_->fsize()),
     c_(c){ assert(c.size()==N_+1); }
  
  Fopr_Chebyshev(const XML::node& cbnode,const Fopr_Herm* Op)
    :Op_(Op),p_(Op_->fsize()),q_(Op_->fsize()),r_(Op_->fsize())
  {
    XML::read(cbnode, "Npoly",N_,MANDATORY);
    c_.resize(N_+1,0.0);
    c_[N_] = 1.0;
  }
  
  const Field mult(const Field& f) const{
    using namespace FieldExpression;

    q_= c_[N_]*f;
    r_= 0.0;

    Field* p = &p_;
    Field* q = &q_;
    Field* r = &r_;
    
    for(int j=0; j<N_; ++j){
      *p = 2.0*Op_->mult(*q) -(*r) +c_[N_-j-1]*f;
      Field* s = r;
      r = q;
      q = p;
      p = s;
    }
    *q -= Op_->mult(*r);
    return *q;
  }

  double func(double x) const{ return Chebyshev::series(Op_->func(x),c_); }
  size_t fsize()const {return Op_->fsize();}
};
#endif
