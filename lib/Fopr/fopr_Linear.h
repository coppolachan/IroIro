#ifndef FOPR_LINEAR_INCLUDED
#define FOPR_LINEAR_INCLUDED

#include "fopr.h"
#include "Fields/field_expressions.hpp"

class Fopr_Linear: public Fopr_Herm{
private:
  const Fopr_Herm* Op_;
  double a_,b_;
public:
  Fopr_Linear(double a,double b,const Fopr_Herm* Op)
    :a_(a),b_(b),Op_(Op){}
  
  Fopr_Linear(XML::node lnode,const Fopr_Herm* Op):Op_(Op){
    XML::read(lnode,"slope",a_,MANDATORY);
    XML::read(lnode,"intercept",b_,MANDATORY);
  }
  
  double func(double x) const{return a_*Op_->func(x) +b_;}
  
  const Field mult(const Field& f) const{
    using namespace FieldExpression;
    
    Field w = Op_->mult(f);
    w *= a_;
    w += b_*f;
    return w;
  }

  size_t fsize()const {return Op_->fsize();}
};

#endif
