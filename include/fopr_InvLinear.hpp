#ifndef FOPR_INVLINEAR_INCLUDED
#define FOPR_INVLINEAR_INCLUDED

#include "fopr.h"
#include "Solver/solver.hpp"
#include "Fields/field_expressions.hpp"

class Fopr_InvLinear:public Fopr_Herm {
private:
  const Fopr_Herm* Op_; /// ensure to be positive definite operator.
  const Solver* slv_;   /// ensure to contain same operator as Op_.
  double a_,b_; 
public:
  Fopr_InvLinear(double a,double b,const Fopr_Herm* Op,const Solver* slv)
    :a_(a),b_(b),Op_(Op),slv_(slv){}

  Fopr_InvLinear(const XML::node& node,const Fopr_Herm* Op,const Solver* slv)
    :Op_(Op),slv_(slv){
    XML::read(node,"slope",a_,MANDATORY);
    XML::read(node,"intercept",b_,MANDATORY);
  }

  const Field mult(const Field& f)const{
    using namespace FieldExpression;

    Field w(f.size());
    SolverOutput so = slv_->solve(w,f);
    w *= a_;
    w += b_*f;
    return w;
  }
  
  double func(double x)const{ return a_/Op_->func(x) +b_; }
  size_t fsize()const{return Op_->fsize();}
};

#endif
