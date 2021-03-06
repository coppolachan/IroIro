#ifndef FOPR_EXP_INCLUDED
#define FOPR_EXP_INCLUDED

#include "fopr.h"
#include "fopr_Linear.h"
#include "Fields/field_expressions.hpp"
#include "pugi_interface.h"
#include "include/timings.hpp"
#include <omp.h>
#include <math.h>

class Fopr_Exp: public Fopr_Herm{
private:
  const Fopr_Herm* Op_;
  int N_;
  double a_;
  double expnt(double x,int n=1)const{
    return N_== n ? 1.0 : 1.0 +a_/n*x*expnt(x,n+1);
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
  double func(double x)const{return expnt(Op_->func(x));}
  size_t fsize()const {return Op_->fsize();}
};


class Fopr_Lexp: public Fopr_Herm{
  const Fopr_Linear* Op_;
  const Fopr_Herm* OpH_;
  int N_;
public:
  Fopr_Lexp(int N,double a,const Fopr_Herm* Op)
    :N_(N),Op_(new Fopr_Linear(a/N,1.0,Op)),OpH_(Op){}

 Fopr_Lexp(XML::node node,const Fopr_Herm* Op):OpH_(Op){
    XML::read(node,"exp_approx",N_,MANDATORY);
    double a;
    XML::read(node,"coefficient",a,MANDATORY);
    Op_= new Fopr_Linear(a/N_,1.0,Op);
  }
  ~Fopr_Lexp(){if(Op_) delete Op_;}
  
  const Field mult(const Field& f)const{
    Field tmp =f;
    Field tmp2=f;
    for(int i=0;i<N_;i++) tmp = Op_->mult(tmp);


    /*
    // Optimizations
    double fake_a = 1.0/N_;
    double *tadd = tmp.getaddr(0);
    double *t2add = tmp2.getaddr(0);
    int Nvol = CommonPrms::instance()->Nvol();
    long double lexp_timer, tmp_timer, inner_timer = 0;
    FINE_TIMING_START(lexp_timer);
    for(int i=0;i<N_;i++) {
      FINE_TIMING_START(tmp_timer);
      tmp2 = OpH_->mult(tmp);
      FINE_TIMING_END(tmp_timer);
      inner_timer += tmp_timer;
      
      
#pragma omp parallel 
      {
	int nid = omp_get_num_threads();
	int is = omp_get_thread_num()*Nvol/nid;
	int ns = Nvol/nid;
	
	BGWilsonLA_MultAddScalar(tadd,t2add,fake_a,ns);
      }
     
 	
    }
    FINE_TIMING_END(lexp_timer);
    CCIO::cout << " Timing - Lexp     :  "<< lexp_timer << " seconds\n";
    CCIO::cout << " Timing - Lexp In  :  "<< inner_timer << " seconds\n";
    */


    return tmp;
  }

  double func(double x)const{return pow(Op_->func(x),N_);}
  size_t fsize()const {return Op_->fsize();}
};

#endif
