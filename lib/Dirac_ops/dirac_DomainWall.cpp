/*!--------------------------------------------------------------------------
 * @file dirac_DomainWall.cpp
 *
 * @brief Definition of class methods for Dirac_optimalDomainWall (5d operator)
 *
 *-------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>

#include "dirac_DomainWall.hpp"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"

using namespace std;

// Constructors for DomainWall Parameters classes
//======================================================================
Dirac_optimalDomainWall_params::
Dirac_optimalDomainWall_params(XML::node DWF_node){
  std::string Precond_string;
  XML::read(DWF_node, "Preconditioning", Precond_string, MANDATORY);
  XML::read(DWF_node, "N5d", N5_, MANDATORY);
  XML::read(DWF_node, "wilson_mass", M0_, MANDATORY);
  XML::read(DWF_node, "b", b_, MANDATORY);
  XML::read(DWF_node, "c", c_, MANDATORY);
  XML::read(DWF_node, "mass", mq_, MANDATORY);

  CCIO::cout<<"mass="<<mq_<<std::endl;
  CCIO::cout<<"Dirac_optimalDomainWall_params::N5_="<<N5_<<std::endl;

  XML::node ApproxNode = DWF_node.child("approximation");
  if(ApproxNode !=NULL) {
    const char* Approx_name = ApproxNode.attribute("name").value();
    if(!strcmp(Approx_name, "Zolotarev")){
      double lambda_min, lambda_max;
      XML::read(ApproxNode, "lambda_min", lambda_min); 
      XML::read(ApproxNode, "lambda_max", lambda_max); 
      omega_= DomainWallFermions::getOmega(N5_,lambda_min,lambda_max);
    }
    if (!strcmp(Approx_name, "Tanh"))  
      for (int s=0; s<N5_; ++s) omega_.push_back(1.0);
  }else{
    CCIO::cout << "Error: missing [approximation] node or wrong entry\n";
    abort();
  }
  // setup of the member arrays
  set_arrays();
  if (!EnumString<Preconditioners>::To( Preconditioning_, Precond_string )){
    CCIO::cerr << "Error: string ["<< Precond_string <<"] not valid" 
	       <<std::endl;
    abort();
  } else {
    CCIO::cout << "Choosing preconditioner type: "
	       << Precond_string << " Code: "<< Preconditioning_ <<std::endl;
  }
}

Dirac_optimalDomainWall_params::
Dirac_optimalDomainWall_params(double b,double c,double M0,double mq,
			       const std::vector<double>& omega, 
			       Preconditioners Preconditioning)
  :N5_(omega.size()),b_(b),c_(c),M0_(M0),mq_(mq),omega_(omega),
   Preconditioning_(Preconditioning){
  set_arrays();
}

void Dirac_optimalDomainWall_params::set_arrays(){

  for(int s=0; s<N5_; ++s){
    bs_.push_back(0.5*(b_*omega_[s] +c_));
    cs_.push_back(0.5*(b_*omega_[s] -c_));
    dp_.push_back(bs_[s]*(4.0 +M0_)+1.0);
    dm_.push_back(1.0 -cs_[s]*(4.0 +M0_));
  }
  es_.push_back(dm_[N5_-1]/dp_[0]);
  fs_.push_back(dm_[0]);
  for (int s=1; s<N5_-1; ++s) {
    es_.push_back(es_[s-1]*(dm_[s-1]/dp_[s]));
    fs_.push_back(fs_[s-1]*(dm_[s]/dp_[s-1]));
  }
}

///////////////////////////////////////////////////////////////////////////
/*! @brief mult in heavy quark limit wihch is useful in preconditionings */

const Field Dirac_optimalDomainWall::mult_hq(const Field& f5) const{
  using namespace FieldExpression;

  assert(f5.size()==fsize_);
  Field w5(fsize_);

  for(int s=0; s<N5_; ++s) {
    Field v = get4d(f5,s);
    v *= Params.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s) {
    Field v = proj_m(get4d(f5,s+1));
    v *= -Params.dm_[s];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s) {
    Field v = proj_p(get4d(f5,s-1));
    v *= -Params.dm_[s];
    add5d(w5,v,s);
  }

  Field v = proj_p(get4d(f5,N5_-1));
  v *= mq_*Params.dm_[0];
  add5d(w5,v,0);

  v = proj_m(get4d(f5,0));
  v *= mq_*Params.dm_[N5_-1];
  add5d(w5,v,N5_-1);

  return w5;
}

const Field Dirac_optimalDomainWall::mult_hq_dag(const Field& f5) const{
  return R5(mult_hq(R5(f5)));
  /*
  using namespace FieldExpression;

  assert(f5.size()==fsize_);
  Field w5(fsize_);

  for(int s=0; s<N5_; ++s){
    Field v = get4d(f5,s);
    v *= Params.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s){
    Field v = proj_p(get4d(f5,s+1));
    v *= -Params.dm_[s+1];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s){
    Field v = proj_m(get4d(f5,s-1));
    v *= -Params.dm_[s-1];
    add5d(w5,v,s);
  }
  Field v = proj_m(get4d(f5,N5_-1));
  v *= mq_*Params.dm_[N5_-1];
  add5d(w5,v,0);

  v = proj_p(get4d(f5,0));
  v *= mq_*Params.dm_[0];
  add5d(w5,v,N5_-1);

  return w5;
  */
}

const Field Dirac_optimalDomainWall::mult_hq_inv(const Field& f5) const{
 using namespace FieldExpression;

  assert(f5.size()==fsize_);
  Field w5(f5);

  for (int s=1; s<N5_; ++s) {
    Field lpf =  proj_p(get4d(w5,s-1));
    lpf *= (Params.dm_[s]/Params.dp_[s-1]);
    add5d(w5,lpf,s);
  }
  for (int s=0; s< N5_-1; ++s) {
    Field ey = proj_m(get4d(w5,s));
    ey *= -mq_* Params.es_[s];
    add5d(w5,ey,N5_-1);
  }

  Field v =  get4d(w5,N5_-1);
  v *= 1.0/(Params.dp_[N5_-1]
            +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]);
  set5d(w5,v,N5_-1);
  for(int s=N5_-2; s>=0; --s) {
    Field lmf = proj_m(get4d(w5,s+1));
    lmf *= Params.dm_[s];
    add5d(w5,lmf,s);
    Field fy = proj_p(get4d(w5,N5_-1));
    fy *= -mq_*Params.fs_[s];
    add5d(w5,fy,s);
    v = get4d(w5,s);
    v*= 1.0/ Params.dp_[s];
    set5d(w5,v,s);
  }
  return w5;
}

const Field Dirac_optimalDomainWall::mult_hq_dinv(const Field& f5) const{
  return R5(mult_hq_inv(R5(f5)));
  /*
  assert(f5.size()==fsize_);   
  Field t5(fsize_);
  t5 = f5;
  
  // LU preconditioning : ((LU)^T)^-1 = (U^T L^T)^-1 = (L^T)^-1 (U^T)^-1
  Field v = get4d(t5,0);
  v *= 1.0/Params.dp_[0];
  set5d(t5,v,0);

  for(int s=1; s<N5_-1; ++s){
    Field lmf = proj_m(get4d(t5,s-1));                                         
    lmf *= Params.dm_[s-1];                                                     
    add5d(t5,lmf,s);                                                            
    v = get4d(t5,s);                                                            
    v*= 1.0/Params.dp_[s];                                                      
    set5d(t5,v,s);                                                              
  }                                                                             
  v = proj_m(get4d(t5,N5_-2));                                                  
  v*= Params.dm_[N5_-2];                                                        
  add5d(t5,v,N5_-1);                                                            
  for(int s=0; s<N5_-1; ++s) {                                                  
    Field fy = proj_p(get4d(t5,s));                                             
    fy *= -mq_*Params.fs_[s];                                            
    add5d(t5,fy,N5_-1);                                                         
  }                                                                             
  v = get4d(t5,N5_-1);                                                          
  v*= 1.0/(Params.dp_[N5_-1]                                                    
           +mq_*Params.dm_[N5_-2]*Params.es_[N5_-2]);                    
  set5d(t5,v,N5_-1);                                                            
                                                                                
  for(int s=N5_-2; s>=0; --s){                                                  
    Field lpf = proj_p(get4d(t5,s+1));                                          
    lpf *= (Params.dm_[s+1]/Params.dp_[s]);                                     
    add5d(t5,lpf,s);                                                            
    Field ey = proj_m(get4d(t5,N5_-1));                                         
    ey *= -mq_*Params.es_[s];                                            
    add5d(t5,ey,s);                                                             
  }                                                                             
  return t5;                                                                    
  */
}

/*! @brief definitions of D_dwf */
void Dirac_optimalDomainWall::mult_a0(Field& w5, const Field& f5) const{ 
  using namespace FieldExpression;
  assert(w5.size()==f5.size());

  for(int s=0; s<N5_; ++s) {
    Field lpf = proj_p(get4d(f5,(s+N5_-1)%N5_));
    if(s==0)     lpf *= -mq_;
    Field lmf = proj_m(get4d(f5,(s+1)%N5_));
    if(s==N5_-1) lmf *= -mq_;

    Field v = get4d(f5,s);
    v *= Params.bs_[s];
    v += Params.cs_[s]*(lpf +lmf);

    Field w = Dw_.mult(v);
    w *= 4.0+M0_;
    set5d(w5,w,s);
  }
}

void Dirac_optimalDomainWall::mult_a1(Field& w5, const Field& f5) const{ 
  using namespace FieldExpression;
  assert(w5.size()==f5.size());
  
  for(int s=0; s<N5_; ++s) {
    Field lpf = proj_p(get4d(f5,(s+N5_-1)%N5_));
    if(s==0)     lpf *= -mq_;
    Field lmf = proj_m(get4d(f5,(s+1)%N5_));
    if(s==N5_-1) lmf *= -mq_;

    Field v = get4d(f5,s);
    v *= Params.bs_[s];
    v += Params.cs_[s]*(lpf +lmf);

    Field w = Dw_.mult(v);          
    w *= 4.0+M0_;
    w += get4d(f5,s) -lpf -lmf;
    set5d(w5,w,s);
  }
}

const Field Dirac_optimalDomainWall::mult(const Field& f5) const{
  Field w5(fsize_);
  (this->*mult_core)(w5,f5);
  return w5;
}

void Dirac_optimalDomainWall::mult_dag_a0(Field& w5,const Field& f5) const{
  using namespace FieldExpression;
  assert(w5.size()==f5.size());
  
  for(int s=0; s<N5_; ++s) {
    Field lpf = proj_p(get4d(f5,(s+N5_-1)%N5_));
    if(s==N5_-1)     lpf *= -mq_;
    Field lmf = proj_m(get4d(f5,(s+1)%N5_));
    if(s==0) lmf *= -mq_;
    Field v = get4d(f5,s);

    v *= Params.bs_[s];
    v += Params.cs_[s]*(lpf +lmf);

    Field w = Dw_.mult_dag(v);          
    w *= 4.0+M0_;
    set5d(w5,w,s);
  }
}

void Dirac_optimalDomainWall::mult_dag_a1(Field& w5,const Field& f5) const{
  using namespace FieldExpression;
  assert(f5.size()==fsize_);
  Field v5(fsize_);
  Field t5(fsize_);
  t5 = f5;
   
  for(int s=0; s<N5_; ++s){
    Field dv = Dw_.mult_dag(get4d(t5,s));
    dv *= (4.0+M0_)*Params.bs_[s];
    set5d(w5,dv,s);
    dv *= (Params.cs_[s]/Params.bs_[s]);
    set5d(v5,dv,s);
  }
  w5 += t5;
  v5 -= t5;
 
  for(int s = 0; s < N5_; ++s){
    Field lpf = proj_p(get4d(v5,(s+1)%N5_));
    if(s == N5_-1) lpf *= -mq_;
    Field lmf = proj_m(get4d(v5,(s+N5_-1)%N5_));
    if(s == 0)     lmf *= -mq_;
    add5d(w5,lpf,s);
    add5d(w5,lmf,s);
  }
}

/////// related useful functions  //////
const Field Dirac_optimalDomainWall::mult_dag(const Field& f5) const{
  //  return R5g5(mult(R5g5(f5)));
  Field w5(fsize_);
  (this->*mult_dag_core)(w5,f5);
  return w5;
}

const Field Dirac_optimalDomainWall::Dminus(const Field& f5) const{
  //1-c_s * D_w(-M)
  Field w5(fsize_);
  w5 = f5;
  for(int s = 0; s < N5_; ++s) {
    Field temp =  Dw_.mult(get4d(f5,s));
    temp *= -Params.cs_[s]; // = [-c_s * D_w(-M)]f5
    add5d(w5, temp, s); //= [1-c_s * D_w(-M)]f5
  }
  return w5;
}

const Field Dirac_optimalDomainWall::gamma5(const Field& f5) const{
  Field w5(fsize_); 
  for(int s = 0; s < N5_; ++s) set5d(w5,gamma5_4d(get4d(f5,s)),s);
  return w5; 
}

const Field Dirac_optimalDomainWall::R5(const Field& f5) const{
  Field w5(fsize_); 
  for(int s = 0; s < N5_; ++s) set5d(w5,get4d(f5,s),N5_-s-1);
  return w5; 
}

const Field Dirac_optimalDomainWall::R5g5(const Field& f5) const{
  //return R5(gamma5(f5)); 
  Field w5(fsize_);
  for(int s=0; s<N5_; ++s) set5d(w5,gamma5_4d(get4d(f5,s)),N5_-s-1);
  return w5;
}

const Field Dirac_optimalDomainWall::Bproj( const Field& f5) const{ 
  Field f4 = proj_p(get4d(f5,N5_-1));
  f4 += proj_m(get4d(f5,0));
  //  Field f4 = proj_p(get4d(f5,0));
  //  f4 += proj_m(get4d(f5,N5_-1));
  return f4;
}

const Field Dirac_optimalDomainWall::Bproj_dag(const Field& f4) const{
  Field f5(fsize_);
  set5d(f5,proj_p(f4),N5_-1);
  set5d(f5,proj_m(f4),0);
  //  set5d(f5,proj_p(f4),0);
  //  set5d(f5,proj_m(f4),N5_-1);
  return f5;
}

const Field Dirac_optimalDomainWall::proj_p(const Field& f4) const{
  Field w4 = Dw_.proj_p(f4);
  return w4;
}

const Field Dirac_optimalDomainWall::proj_m(const Field& f4) const{
  Field w4 = Dw_.proj_m(f4);
  return w4;
}

//////////////////////////////////////////////////////////////////
/*! @brief contribution to the MD-force from forward difference */
void Dirac_optimalDomainWall::
md_force_p(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;

  for(int s=0; s<N5_; ++s){
    Field lpf = proj_p(get4d(phi,(s+N5_-1)%N5_));
    if(s == 0)     lpf *= -mq_;
    Field lmf = proj_m(get4d(phi,(s+1    )%N5_));
    if(s == N5_-1) lmf *= -mq_;

    Field w = get4d(phi,s); 

    w *= Params.bs_[s];
    w += Params.cs_[s]*(lpf +lmf);

    Dw_.md_force_p(fce,w,get4d(psi,s));
  }
}  

/*! @brief contribution to the MD-force from backward difference */
void Dirac_optimalDomainWall::
md_force_m(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;

  for(int s=0; s<N5_; ++s){
    Field lpf = proj_p(get4d(phi,(s+N5_-1)%N5_));
    if(s == 0)     lpf *= -mq_;
    Field lmf = proj_m(get4d(phi,(s+1    )%N5_));
    if(s == N5_-1) lmf *= -mq_;

    Field w = get4d(phi,s); 

    w *= Params.bs_[s];
    w += Params.cs_[s]*(lpf +lmf);

    Dw_.md_force_m(fce,w,get4d(psi,s));
  }
}  

/*! @brief total MD-force */
const Field Dirac_optimalDomainWall::
md_force(const Field& phi,const Field& psi) const{
  Field fce(gsize_);
  md_force_p(fce,phi,psi);
  md_force_m(fce,phi,psi);
  fce *= -0.5;
  return fce;
}
//
Preconditioner* Dirac_optimalDomainWall::
choose_Preconditioner(int PrecondID){
  switch (PrecondID){
  case NoPreconditioner:
    return new NoPrecond(this);
  case LUPreconditioner:
    return new LUPrecond(this);
  default:
    return new NoPrecond(this);
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
  @brief Namespace definining useful functions for DomainWallFermions
*/
namespace DomainWallFermions {

  inline double set_vs( int is, int ns, double kprime ){
    double ekprime = gsl_sf_ellint_Kcomp( kprime , 0 ); 
    double vs = is * ekprime / ns; 
    return vs;
  }

  double read_wilson_mass(const XML::node& node){
    double M0;
    XML::read(node,"wilson_mass", M0);
    return M0;
  }

  const vector<double> getOmega(int Ns,double lambda_min,double lambda_max){
    double u, m;
    double sn, cn, dn; 
    double kprime = sqrt(1.0-(lambda_min/lambda_max)*(lambda_min/lambda_max));
    
    vector<double> omegas(Ns);

    for(int ii=0; ii<Ns; ++ii){
      int is = 2*ii + 1;
      m = kprime*kprime;
      double vs = set_vs( is, Ns*2, kprime );
      gsl_sf_elljac_e( vs, m, &sn, &cn, &dn );
      double sn2 = sn * sn;
      double kappaprime2 = kprime * kprime;
      omegas[ii] = (1.0/lambda_min)* sqrt(1.0-kappaprime2*sn2);
    }
#if VERBOSITY>2
    for( int ii=0; ii<Ns; ++ii) printf("%24.16E\n", omegas[ii] );
#endif
    return omegas;
  }
}
