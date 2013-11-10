/*!--------------------------------------------------------------------------
 * @file dirac_DomainWall.cpp
 * @brief Definition of class methods for Dirac_optimalDomainWall (5d operator)
 Time-stamp: <2013-10-17 15:38:46 cossu>
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
Dirac_optimalDomainWall_params(XML::node DWF_node,DWFType Type){
  /* temporal hack*/
  XML::node mynode = DWF_node;
  XML::descend(mynode,"BaseKernel", MANDATORY);
  XML::read(mynode, "mass", M0_, MANDATORY);

  XML::read(DWF_node, "N5d", N5_, MANDATORY);
  XML::read(DWF_node, "b", b_, MANDATORY);
  XML::read(DWF_node, "c", c_, MANDATORY);
  XML::read(DWF_node, "mass", mq_, MANDATORY);

  if(Type == PauliVillars) mq_= 1.0;
  
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
      for(int s=0; s<N5_; ++s) omega_.push_back(1.0);
  }else{
    CCIO::cout << "Error: missing [approximation] node or wrong entry\n";
    abort();
  }
  set_arrays();   // setup of the member arrays
}

Dirac_optimalDomainWall_params::
Dirac_optimalDomainWall_params(double b,double c,double M0,double mq,
			       const std::vector<double>& omega)
  :N5_(omega.size()),b_(b),c_(c),M0_(M0),mq_(mq),omega_(omega){
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


/* @brief Namespace definining useful functions for DomainWallFermions*/
namespace DomainWallFermions {

  inline double set_vs( int is, int ns, double kprime ){
    double ekprime = gsl_sf_ellint_Kcomp( kprime , 0 ); 
    double vs = is * ekprime / ns; 
    return vs;
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

const Field Dirac_optimalDomainWall::mult(const Field& f5) const{
  Field w5(f5size_);
  (this->*mult_core)(w5,f5);
  return w5;
}

const Field Dirac_optimalDomainWall::mult_dag(const Field& f5) const{
  Field w5(f5size_);
  (this->*mult_dag_core)(w5,f5);
  return w5;
}

/*! @brief total MD-force */
const Field Dirac_optimalDomainWall::
md_force(const Field& phi,const Field& psi) const{

  Field fce(Dw_->gsize());
  md_force_p(fce,phi,psi);
  md_force_m(fce,phi,psi);
  fce *= -0.5;
  return fce;
}

const Field Dirac_optimalDomainWall::Dminus(const Field& f5) const{
  //1-c_s * D_w(-M)
  Field w5(f5size_);
  w5 = f5;
  for(int s=0; s<N5_; ++s) {
    Field temp =  Dw_->mult(get4d(f5,s));
    temp *= -Params_.cs_[s]; // = [-c_s * D_w(-M)]f5
    add5d(w5,temp,s); //= [1-c_s * D_w(-M)]f5
  }
  return w5;
}

const Field Dirac_optimalDomainWall::gamma5(const Field& f5) const{
  Field w5(f5size_); 
  
  for(int s=0; s<N5_; ++s){
    for(int site=0; site<Nvol_; ++site){
      gamma5core(w5.getaddr(ff_.index(0,site,s)),
		 const_cast<Field&>(f5).getaddr(ff_.index(0,site,s)));
    }
  }
  return w5; 
}

const Field Dirac_optimalDomainWall::R5(const Field& f5) const{
  Field w5(f5size_); 
  for(int s=0; s<N5_; ++s) set5d(w5,get4d(f5,s),N5_-s-1);
  return w5; 
}

const Field Dirac_optimalDomainWall::R5g5(const Field& f5) const{
  Field w5(f5size_);

  for(int s=0; s<N5_; ++s){
    for(int site=0; site<Nvol_; ++site){
      gamma5core(w5.getaddr(ff_.index(0,site,N5_-s-1)),
		 const_cast<Field&>(f5).getaddr(ff_.index(0,site,s)));
    }
  }
  return w5;
}
/*
const Field Dirac_optimalDomainWall::Bproj( const Field& f5) const{ 
  Field f4(f4size_),t4(f4size_);

  proj_p(f4,f5,N5_-1);
  proj_m(t4,f5,0);
  f4 += t4;
  return f4;
}

const Field Dirac_optimalDomainWall::Bproj_dag(const Field& f4) const{
  Field f5(f5size_),t4(f4size_);
  proj_p(t4,f4);
  set5d(f5,t4,N5_-1);
  proj_m(t4,f4);
  set5d(f5,t4,0);
  //  t4 -= f4;
  //  set5d_c(f5,t4,-1.0,0);
  return f5;
}
*/
void Dirac_optimalDomainWall::proj_p(Field& w,const Field& f5,int s)const{
  for(int site=0; site<Nvol_; ++site)
    projPcore(w.getaddr(ff_.index(0,site)),
	      const_cast<Field&>(f5).getaddr(ff_.index(0,site,s)));
}

void Dirac_optimalDomainWall::proj_m(Field& w,const Field& f5,int s)const{
  for(int site=0; site<Nvol_; ++site)
    projMcore(w.getaddr(ff_.index(0,site)),
	      const_cast<Field&>(f5).getaddr(ff_.index(0,site,s)));
}

void Dirac_optimalDomainWall::gamma5_4d(Field& w,const Field& f)const{
  w = Dw_->gamma5(f);
}

const Field Dirac_optimalDomainWall::get4d(const Field& f5,int s) const{
  Field w(f4size_);
  for(int i=0; i<f4size_; ++i) w.set(i, f5[s*f4size_+i]);
  return w;
}
void Dirac_optimalDomainWall::add5d(Field& f5,const Field& f4_1,const Field& f4_2,int s) const{
  for(int i=0; i<f4size_; ++i) f5.add(s*f4size_+i,f4_1[i]+f4_2[i]);
}
void Dirac_optimalDomainWall::add5d(Field& f5,const Field& f4,int s) const{
  f5.add(std::slice(s*f4size_,f4size_,1),f4.getva());
}
void Dirac_optimalDomainWall::set5d(Field& f5,const Field& f4,int s) const{
  for(int i=0; i<f4size_; ++i) f5.set(s*f4size_+i,f4[i]);
}
void Dirac_optimalDomainWall::mul5d(Field& f5,double fac,int s) const{
  for(int i=0; i<f4size_; ++i) f5.set(s*f4size_+i,fac*f5[s*f4size_+i]);
}


//////////// environment-dependent parts /////////////////////////
#ifdef IBM_BGQ_WILSON
#include "dirac_DomainWall_BGQ.code"
#include "domainWallSolver_BGQ.code"
#else 
/*! @brief mult without 4D-hopping parameters */

const Field Dirac_optimalDomainWall::mult_hop5(const Field& f5) const{
  Field w5(f5size_), v(f4size_);

  for(int s=0; s<N5_; ++s){
    Field v = get4d(f5,s);
    v *= Params_.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s){
    proj_m(v,f5,s+1);
    v *= -Params_.dm_[s];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s){
    proj_p(v,f5,s-1);
    v *= -Params_.dm_[s];
    add5d(w5,v,s);
  }

  proj_p(v,f5,N5_-1);
  v *= mq_*Params_.dm_[0];
  add5d(w5,v,0);

  proj_m(v,f5,0);
  v *= mq_*Params_.dm_[N5_-1];
  add5d(w5,v,N5_-1);
  return w5;
}

const Field Dirac_optimalDomainWall::mult_hop5_dag(const Field& f5) const{
  Field w5(f5size_), v(f4size_);

  assert(f5.size()==f5size_);
  for(int s=0; s<N5_; ++s){
    v = get4d(f5,s);
    v *= Params_.dp_[s];
    set5d(w5,v,s);
  }
  for(int s=0; s<N5_-1; ++s){
    proj_p(v,f5,s+1);
    v *= -Params_.dm_[s+1];
    add5d(w5,v,s);
  }
  for(int s=1; s<N5_; ++s){
    proj_m(v,f5,s-1);
    v *= -Params_.dm_[s-1];
    add5d(w5,v,s);
  }
  proj_m(v,f5,N5_-1);
  v *= mq_*Params_.dm_[N5_-1];
  add5d(w5,v,0);

  proj_p(v,f5,0);
  v *= mq_*Params_.dm_[0];
  add5d(w5,v,N5_-1);
  return w5;
}

const Field Dirac_optimalDomainWall::mult_hop5_inv(const Field& f5) const{
  Field w5(f5);
  Field lf(f4size_),fy(f4size_);

  proj_m(fy,w5,0);
  fy *= -mq_* Params_.es_[0];
  add5d(w5,fy,N5_-1);

  for(int s=1; s<N5_-1; ++s){
    proj_p(lf,w5,s-1);
    proj_m(fy,w5,s);

    lf *= Params_.dm_[s]/Params_.dp_[s-1];
    add5d(w5,lf,s);

    fy *= -mq_*Params_.es_[s];
    add5d(w5,fy,N5_-1);
  }
  proj_p(lf,w5,N5_-2);
  lf *= (Params_.dm_[N5_-1]/Params_.dp_[N5_-2]);
  add5d(w5,lf,N5_-1);

  mul5d(w5,1.0/(Params_.dp_[N5_-1] +mq_*Params_.dm_[N5_-2]*Params_.es_[N5_-2]),
	N5_-1);  

  for(int s=N5_-2; s>=0; --s){
    proj_m(lf,w5,s+1);
    lf *= Params_.dm_[s];
    add5d(w5,lf,s);
    proj_p(fy,w5,N5_-1);
    fy *= -mq_*Params_.fs_[s];
    add5d(w5,fy,s);
    Field v = get4d(w5,s);
    v*= 1.0/ Params_.dp_[s];
    set5d(w5,v,s);
  }
  return w5;
}

const Field Dirac_optimalDomainWall::mult_hop5_dinv(const Field& f5) const{

  Field w5(f5),lpf(f4size_),ey(f4size_),lmf(f4size_),fy(f4size_),v(f4size_);
  double* lpf_ptr = lpf.getaddr(0);
  double* ey_ptr = ey.getaddr(0);

  v = get4d(w5,0);
  v *= 1.0/Params_.dp_[0];
  set5d(w5,v,0);

  for(int s=1; s<N5_-1; ++s){
    proj_m(lmf,w5,s-1);                                         
    lmf *= Params_.dm_[s-1];                                                    

    add5d(w5,lmf,s);                                                            
    v = get4d(w5,s);                                                            
    v*= 1.0/Params_.dp_[s];                                                    
 
    set5d(w5,v,s);                                                              
  }                                                                             
  proj_m(v,w5,N5_-2);                                                  
  v*= Params_.dm_[N5_-2];                                                      
 
  add5d(w5,v,N5_-1);                                                            
  for(int s=0; s<N5_-1; ++s) {                                                  
    proj_p(fy,w5,s);                                             
    fy *= -mq_*Params_.fs_[s];                                            
    add5d(w5,fy,N5_-1);                                                         
  }                                                                             
  v = get4d(w5,N5_-1);                                                          
  v*= 1.0/(Params_.dp_[N5_-1] +mq_*Params_.dm_[N5_-2]*Params_.es_[N5_-2]);    
  set5d(w5,v,N5_-1);                                                            
                                                                                
  for(int s=N5_-2; s>=0; --s){ 
    double* w5_ptr   = w5.getaddr(s*f4size_);
    double fact_lpf = (Params_.dm_[s+1]/Params_.dp_[s]);
    double fact_ey =  mq_*Params_.es_[s];

    proj_p(lpf,w5,  s+1);                                          
    proj_m( ey,w5,N5_-1); 

    for(int i=0; i<f4size_; ++i)
      w5_ptr[i] += fact_lpf*lpf_ptr[i]- fact_ey*ey_ptr[i];
  }
  return w5;    
}

/*! @brief definitions of D_dwf */
void Dirac_optimalDomainWall::mult_full(Field& w5, const Field& f5) const{ 

  Field v(f4size_),lpf(f4size_), lmf(f4size_),w(f4size_);
  double* v_ptr   = v.getaddr(0);
  double* lpf_ptr = lpf.getaddr(0);
  double* lmf_ptr = lmf.getaddr(0);
  double* w_ptr = w.getaddr(0);
  double mass_fact= 4.0+M0_;

  for(int s=0; s<N5_; ++s){
    double* f5_ptr = const_cast<Field&>(f5).getaddr(s*f4size_);
    double* w5_ptr = w5.getaddr(s*f4size_);

    proj_p(lpf,f5,(s+N5_-1)%N5_);
    if(s==0)     lpf *= -mq_;
    proj_m(lmf,f5,(s+1)%N5_);
    if(s==N5_-1) lmf *= -mq_;

    for(int i=0; i<f4size_; ++i){
      lpf_ptr[i] += lmf_ptr[i];
      v_ptr[i] = Params_.bs_[s]*f5_ptr[i]+Params_.cs_[s]*lpf_ptr[i];
    }
    w = Dw_->mult(v);
    for(int i=0; i<f4size_; ++i)
      w5_ptr[i] = mass_fact*w_ptr[i]+ f5_ptr[i] -lpf_ptr[i];
  }
}

void Dirac_optimalDomainWall::mult_dag_full(Field& w5,const Field& f5) const{
  assert(w5.size()==f5.size());
  Field v5(f5size_);
  Field lpf(f4size_),lmf(f4size_),w(f4size_);

  for(int s=0; s< N5_; ++s){
    int spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    double bs = (4.0+M0_)*Params_.bs_[s];
    double cs = (4.0+M0_)*Params_.cs_[s];

    Field w = Dw_->mult_dag(get4d(f5,s));
    double* w_ptr = w.getaddr(0);
    
    for(int i=0; i<f4size_; ++i){
      w5_ptr[i] = bs*w_ptr[i];
      v5_ptr[i] = cs*w_ptr[i];
    }
    // do not change this
    for(int i=0; i<f4size_; ++i){
      w5_ptr[i] += f5_ptr[i];
      v5_ptr[i] -= f5_ptr[i];
    }
  }
  for(int s = 0; s < N5_; ++s){
    proj_p(lpf,v5,(s+1)%N5_);
    if(s == N5_-1) lpf *= -mq_;
    proj_m(lmf,v5,(s+N5_-1)%N5_);
    if(s == 0)     lmf *= -mq_;
    add5d(w5,lpf,lmf,s);
  }
}

void Dirac_optimalDomainWall::mult_offdiag(Field& w5,const Field& f5) const{ 
  using namespace FieldExpression;
  Field lpf(f4size_),lmf(f4size_),w(f4size_),v(f4size_);
  double mass_fact= 4.0+M0_;

  for(int s=0; s<N5_; ++s) {
    proj_p(lpf,f5,(s+N5_-1)%N5_);
    if(s==0)     lpf *= -mq_;
    proj_m(lmf,f5,(s+1)%N5_);
    if(s==N5_-1) lmf *= -mq_;

    lpf += lmf;
    v = get4d(f5,s);
    v *= Params_.bs_[s];
    v += Params_.cs_[s]*lpf;
    w = Dw_->mult(v);
    w *= mass_fact;
    set5d(w5,w,s);
  }
}

void Dirac_optimalDomainWall::mult_dag_offdiag(Field& w5,const Field& f5) const{
  assert(w5.size()==f5.size());
  Field v5(f5size_),lpf(f4size_),lmf(f4size_);

  for(int s=0; s<N5_; ++s){
    int spin_idx = s*f4size_;
    double* f5_ptr = const_cast<Field&>(f5).getaddr(spin_idx);
    double* w5_ptr = w5.getaddr(spin_idx);
    double* v5_ptr = v5.getaddr(spin_idx);
    
    double bs = (4.0+M0_)*Params_.bs_[s];
    double cs = (4.0+M0_)*Params_.cs_[s];
    
    Field w = Dw_->mult_dag(get4d(f5,s));
    double* w_ptr = w.getaddr(0);
    
    for (int i=0; i<f4size_; i++){
      w5_ptr[i] = bs*w_ptr[i];
      v5_ptr[i] = cs*w_ptr[i];
    }
  }
  for(int s=0; s<N5_; ++s){
    proj_p(lpf,v5,(s+1)%N5_);
    if(s == N5_-1) lpf *= -Params_.mq_;
    proj_m(lmf,v5,(s+N5_-1)%N5_);
    if(s == 0)     lmf *= -Params_.mq_;
    add5d(w5,lpf,lmf,s);
  }
}

/*! @brief contribution to the MD-force from forward difference */
void Dirac_optimalDomainWall::
md_force_p(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;
  CCIO::cout << "Calling md_force_p\n";
  Field lpf(f4size_), lmf(f4size_);
  Field xie(f4size_);
  Field w(f4size_);
  double* lpf_ptr  = lpf.getaddr(0);
  double* lmf_ptr  = lmf.getaddr(0);
  double* xie_ptr  = xie.getaddr(0);
  
  for(int s=0; s<N5_; ++s){
    proj_p(lpf,phi,(s+N5_-1)%N5_);
    if(s == 0)     lpf *= -mq_;
    proj_m(lmf,phi,(s+1    )%N5_);
    if(s == N5_-1) lmf *= -mq_;

    Field w = get4d(phi,s); 
    
    w *= Params_.bs_[s];
    w += Params_.cs_[s]*(lpf +lmf);
    
    Dw_->md_force_p(fce,w,get4d(psi,s));  
    
    
  }
}  

/*! @brief contribution to the MD-force from backward difference */
void Dirac_optimalDomainWall::
md_force_m(Field& fce,const Field& phi,const Field& psi)const{
  using namespace FieldExpression;
  Field lpf(f4size_), lmf(f4size_);

  for(int s=0; s<N5_; ++s){
    proj_p(lpf,phi,(s+N5_-1)%N5_);
    if(s == 0)     lpf *= -mq_;
    proj_m(lmf,phi,(s+1    )%N5_);
    if(s == N5_-1) lmf *= -mq_;

    Field w = get4d(phi,s); 

    w *= Params_.bs_[s];
    w += Params_.cs_[s]*(lpf +lmf);

    Dw_->md_force_m(fce,w,get4d(psi,s));
  }
}  

#endif






