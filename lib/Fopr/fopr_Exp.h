#ifndef FOPR_EXP_INCLUDED
#define FOPR_EXP_INCLUDED

#include "fopr.h"
#include "fopr_Linear.h"
#include "Fields/field_expressions.hpp"
#include "pugi_interface.h"
#include <math.h>

class Fopr_Exp: public Fopr_Herm{
private:
  const Fopr_Herm* Op_;
  int N_;
  double a_;
  double expnt(double x,int n=1)const{
    return N_== n ? 1.0 : 1.0 +a_/n*Op_->func(x)*expnt(x,n+1);
  }
  const Field expnt(const Field& f,int n=1)const{
    using namespace FieldExpression;
    Field ff(f);
    if(n!=N_) ff += a_/n*Op_->mult(expnt(f,n+1));
    return ff;
  }
public:
  Fopr_Exp(int N,double a,const Fopr_Herm* Op)
    :N_(N),a_(a),Op_(Op){}

  Fopr_Exp(XML::node node,const Fopr_Herm* Op):Op_(Op){
    XML::read(node,"exp_approx",N_,MANDATORY);
    XML::read(node,"coefficient",a_,MANDATORY);
    assert(Op_!=NULL);
  }
  const Field mult(const Field& f)const{ return expnt(f);}    
  double func(double x)const{return expnt(x);}
  size_t fsize()const {return Op_->fsize();}
};


class Fopr_Lexp: public Fopr_Herm{
  const Fopr_Linear* Op_;
  int N_;
  const Field prod(const Field& f,int n=0)const{
    return N_== n ? f : Op_->mult(prod(f,n+1));
  }
public:
  Fopr_Lexp(int N,double a,const Fopr_Herm* Op)
    :N_(N),Op_(new Fopr_Linear(a/N,1.0,Op)){}

  Fopr_Lexp(XML::node node,const Fopr_Herm* Op){
    XML::read(node,"exp_approx",N_,MANDATORY);
    double a;
    XML::read(node,"coefficient",a,MANDATORY);
    Op_= new Fopr_Linear(a/N_,1.0,Op);
  }
  ~Fopr_Lexp(){if(Op_) delete Op_;}
  
  const Field mult(const Field& f)const{ return prod(f);}    
  double func(double x)const{return pow(Op_->func(x),N_);}
  size_t fsize()const {return Op_->fsize();}
};

#endif
