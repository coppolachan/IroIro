/*!
 * @file dirac_DomainWall.cpp
 *
 * @brief Definition of class methods for Dirac_optimalDomainWall (5d operator)
 */
#include "dirac_DomainWall.hpp"
#include "Communicator/comm_io.hpp"
#include<stdlib.h>
#include<stdio.h>
#include <cassert>

#include<math.h>

#include<gsl/gsl_sf_ellint.h>
#include<gsl/gsl_sf_elljac.h>

using namespace std;

// Constructors
//-------------------------------------------------------------------------------------
Dirac_optimalDomainWall_params::Dirac_optimalDomainWall_params(XML::node DWF_node){
    std::string Precond_string;
    XML::read(DWF_node, "Preconditioning", Precond_string, MANDATORY);
    XML::read(DWF_node, "N5d", N5dim_, MANDATORY);
    XML::read(DWF_node, "b", b_, MANDATORY);
    XML::read(DWF_node, "c", c_, MANDATORY);
    XML::read(DWF_node, "mass", mq_, MANDATORY);
    omega_.resize(N5dim_);
    XML::node ApproxNode = DWF_node.child("approximation");
    if (ApproxNode !=NULL) {
      const char* Approx_name = ApproxNode.attribute("name").value();
      if (!strcmp(Approx_name, "Zolotarev")){
	double lambda_min, lambda_max;
	XML::read(ApproxNode, "lambda_min", lambda_min); 
	XML::read(ApproxNode, "lambda_max", lambda_max); 
	omega_ = DomainWallFermions::getOmega(N5dim_,lambda_min,lambda_max);
      }
      if (!strcmp(Approx_name, "Tanh")){
	for (int s = 0; s < N5dim_; ++s) omega_[s] = 1.0;
      }
    } else {
      CCIO::cout << "Error: missing [approximation] node or wrong entry" << std::endl;
      abort();
    }

    bs_.resize(N5dim_);
    cs_.resize(N5dim_);
    dp_.resize(N5dim_);
    dm_.resize(N5dim_);

    if (!EnumString<Preconditioners>::To( Preconditioning_, Precond_string )){
      CCIO::cerr << "Error: string ["<< Precond_string <<"] not valid"  << std::endl;
      abort();
    } else {
      CCIO::cout << "Choosing preconditioner type: "<< 
	Precond_string << " Code: "<< Preconditioning_ << std::endl;
    }

  }

Dirac_optimalDomainWall_params::Dirac_optimalDomainWall_params(const double b,
							       const double c,
							       const double mass,
							       const std::vector<double> omega){
  b_ = b;
  c_ = c;
  mq_ = mass;
  omega_ = omega;
  bs_.resize(omega.size());
  cs_.resize(omega.size());
  dp_.resize(omega.size());
  dm_.resize(omega.size());
}

Dirac_optimalDomainWall_params::Dirac_optimalDomainWall_params(const Dirac_optimalDomainWall_params& Par,
							       DWFType Type){
  b_ = Par.b_;
  c_ = Par.c_;
  omega_ = Par.omega_;
  Preconditioning_ = Par.Preconditioning_;
  bs_ = Par.bs_;
  cs_ = Par.cs_;
  dp_ = Par.dp_;
  dm_ = Par.dm_;
  
  switch (Type){
  case Standard:
    mq_ = Par.mq_;
    break;
  case PauliVillars:
    mq_ = 1.0;
    break;
    default:
      abort();
  }
}

//--------------------------------------------------------------------------------------

Dirac_optimalDomainWall::Dirac_optimalDomainWall(XML::node DWF_node,
						 const Dirac_Wilson* Kernel)
  :Params(Dirac_optimalDomainWall_params(DWF_node)),
   Dw_(Kernel),
   N5_(Params.omega_.size()),
   f4size_(Dw_->fsize()),
   fsize_(f4size_*N5_),
   gsize_(Dw_->gsize()),
   M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){
  
  for (int s = 0; s < N5_; ++s) {
    Params.bs_[s] = (Params.b_*Params.omega_[s]+Params.c_)/2.0;
    Params.cs_[s] = (Params.b_*Params.omega_[s]-Params.c_)/2.0;
    Params.dp_[s] = Params.bs_[s]*(4.0+M0_)+1.0;
    Params.dm_[s] = 1.0-Params.cs_[s]*(4.0+M0_);
  }
  choose_Preconditioner(Params.Preconditioning_);
  
}

Dirac_optimalDomainWall::Dirac_optimalDomainWall(Dirac_optimalDomainWall_params Prm,
						 const Dirac_Wilson* Kernel)
  :Params(Prm),
   Dw_(Kernel),
   N5_(Params.omega_.size()),
   f4size_(Dw_->fsize()),
   fsize_(f4size_*N5_),
   gsize_(Dw_->gsize()),
   M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){
  choose_Preconditioner(Params.Preconditioning_);
} 

Dirac_optimalDomainWall::Dirac_optimalDomainWall(const double b,
						 const double c,
						 const double mq,
						 const std::vector<double>& omega,
						 Preconditioners Precond,
						 const Dirac_Wilson* Kernel)
  :Params(b,c,mq,omega),
   Dw_(Kernel),
   N5_(Params.omega_.size()),
   f4size_(Dw_->fsize()),
   fsize_(f4size_*N5_),
   gsize_(Dw_->gsize()),
   M0_(1.0/(2.0*(Dw_->getKappa()))-4.0){

    for (int s = 0; s < omega.size(); ++s) {
      Params.bs_[s] = (b*omega[s]+c)/2.0;
      Params.cs_[s] = (b*omega[s]-c)/2.0;
      Params.dp_[s] = Params.bs_[s]*(4.0+M0_)+1.0;
      Params.dm_[s] = 1.0-Params.cs_[s]*(4.0+M0_);
    }

  choose_Preconditioner(Precond);
}


int Dirac_optimalDomainWall::choose_Preconditioner(int PrecondID) {
  switch (PrecondID) {
  case NoPreconditioner:
    Precond_ = new NoPrecond(this);
    return 0;
    break;
  case LUPreconditioner:
    Precond_ = new LUPrecond(this);
    return 0;
    break;
  default:
    Precond_ = new NoPrecond(this);
    return 1;
    break;
  }
}



const Field Dirac_optimalDomainWall::NoPrecond::mult(const Field& f5) const{
  return DWF_->mult(f5);
}

const Field Dirac_optimalDomainWall::NoPrecond::mult_dag(const Field& f5) const{
  return DWF_->mult_dag(f5);
}

const Field Dirac_optimalDomainWall::LUPrecond::mult(const Field& f5) const{
  Field w5(DWF_->fsize_);
  w5 = DWF_->mult(f5);
  
  for (int s=1; s<DWF_->N5_; ++s) {
    Field lpf = DWF_->proj_p(DWF_->get4d(w5,s-1));
    lpf *= (DWF_->Params.dm_[s]/DWF_->Params.dp_[s-1]);
    DWF_->add5d(w5,lpf,s);
  };
  Field v = DWF_->get4d(w5,DWF_->N5_-1);
  v *= 1.0/DWF_->Params.dp_[DWF_->N5_-1];
  DWF_->set5d(w5,v,DWF_->N5_-1);
  for (int s=DWF_->N5_-2; s>=0; --s) {
    Field lmf = DWF_->proj_m(DWF_->get4d(w5,s+1));
    lmf *= DWF_->Params.dm_[s];
    DWF_->add5d(w5,lmf,s);
    v = DWF_->get4d(w5,s);
    v *= 1.0/DWF_->Params.dp_[s];
    DWF_->set5d(w5,v,s);
  }
  return w5;
}

const Field Dirac_optimalDomainWall::LUPrecond::mult_dag(const Field& f5) const{
  assert(f5.size()==DWF_->fsize_);
  Field t5(DWF_->fsize_);
  t5 = f5;
  
  // LU preconditioning : ((LU)^T)^-1 = (U^T L^T)^-1 = (L^T)^-1 (U^T)^-1
  Field v = DWF_->get4d(t5,0);
  v *= 1.0/DWF_->Params.dp_[0];
  DWF_->set5d(t5,v,0);
  for (int s=1; s<DWF_->N5_; ++s) {
    Field lmf = DWF_->proj_m(DWF_->get4d(t5,s-1));
    lmf *= DWF_->Params.dm_[s-1];
    DWF_->add5d(t5,lmf,s);
    v = DWF_->get4d(t5,s);
    v *= 1.0/DWF_->Params.dp_[s];
    DWF_->set5d(t5,v,s);
  }
  for (int s=DWF_->N5_-2; s>=0; --s) {
    Field lpf = DWF_->proj_p(DWF_->get4d(t5,s+1));
    lpf *= (DWF_->Params.dm_[s+1]/DWF_->Params.dp_[s]);
    DWF_->add5d(t5,lpf,s);
  }
  // end precond LU

  //Multiply D_dwf
  return DWF_->mult_dag(t5); 
}


const Field Dirac_optimalDomainWall::mult(const Field& f5) const{ 
  using namespace FieldExpression;

  assert(f5.size()==fsize_);
  Field w5(fsize_);

  for(int s=0; s<N5_; ++s) {
    Field lpf = proj_p(get4d(f5,(s+N5_-1)%N5_));
    if(s == 0)     lpf *= -Params.mq_;
    Field lmf = proj_m(get4d(f5,(s+1)%N5_));
    if(s == N5_-1) lmf *= -Params.mq_;
      
    Field w = get4d(f5,s);
    w -= lpf;
    w -= lmf;
      
    Field v = get4d(f5,s);
    v += (Params.cs_[s]/Params.bs_[s])*lpf;
    v += (Params.cs_[s]/Params.bs_[s])*lmf;

    w += (4.0+M0_)*Params.bs_[s]*(Dw_->mult(v));      
    set5d(w5,w,s);
    }
  return w5;
}


const Field Dirac_optimalDomainWall::mult_dag(const Field& f5) const{
  assert(f5.size()==fsize_);
  Field v5(fsize_);
  Field w5(fsize_);

  Field t5(fsize_);
  t5 = f5;
  
  for(int s=0; s<N5_; ++s){
    Field dv = Dw_->mult_dag(get4d(t5,s));
    dv *= (4.0+M0_)*Params.bs_[s];
    set5d(w5,dv,s);
    dv *= (Params.cs_[s]/Params.bs_[s]);
    set5d(v5,dv,s);
  }
  w5 += t5;
  v5 -= t5;

  for(int s = 0; s < N5_; ++s){
    Field lpf = proj_p(get4d(v5,(s+1)%N5_));
    if(s == N5_-1) lpf *= -Params.mq_;
    Field lmf = proj_m(get4d(v5,(s+N5_-1)%N5_));
    if(s == 0)     lmf *= -Params.mq_;
    add5d(w5,lpf,s);
    add5d(w5,lmf,s);

  }
  return w5;
}

const Field Dirac_optimalDomainWall::
md_force(const Field& phi,const Field& psi) const{

  using namespace FieldExpression;

  Field w5(fsize_);
  for(int s = 0; s < N5_; ++s){

    Field lpf = proj_p(get4d(phi,(s +1)%N5_));
    if(s == N5_-1) lpf *= -Params.mq_;
    Field lmf = proj_m(get4d(phi,(s +N5_-1)%N5_));
    if(s == 0) lmf *= -Params.mq_;

    Field w4 = Field(get4d(phi,s));
    w4 += Params.c_*lpf;
    w4 += Params.c_*lmf;
    
    set5d(w5,w4,s);
  }

  Field force(gsize_);
  for(int s = 0; s < N5_; ++s)
    force += (4.0+M0_)*Params.omega_[s]*Dw_->md_force(get4d(w5,s),get4d(psi,s));

  return force;
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
  return R5(gamma5(f5)); 
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
  Field w4 = Dw_->proj_p(f4);
  return w4;
}

const Field Dirac_optimalDomainWall::proj_m(const Field& f4) const{
  Field w4 = Dw_->proj_m(f4);
  return w4;
}

namespace DomainWallFermions {

  inline double set_vs( int is , int ns , double kprime )
  {
    double ekprime = gsl_sf_ellint_Kcomp( kprime , 0 ); 
    double vs = is * ekprime / ns; 
    return vs;
  }
  

  vector<double> getOmega(int Ns, 
			  double lambda_min, 
			  double lambda_max)
  {
    double u, m;
    double sn , cn , dn; 
    double kprime = sqrt( 1.0 - (lambda_min/lambda_max) *
			  (lambda_min/lambda_max));
    
    vector<double> omegas(Ns);

    for( int ii = 0 ; ii < Ns ; ii ++ )
      {
	int is = 2 * ii + 1;
	m = kprime * kprime;
	double vs = set_vs( is , Ns*2 , kprime );
	gsl_sf_elljac_e( vs , m , &sn , &cn , &dn );
	double sn2 = sn * sn;
	double kappaprime2 = kprime * kprime;
	omegas[ii] = ( 1.0 / lambda_min ) * sqrt( 1.0 - kappaprime2 * sn2 );
      }
    
    #ifdef VERBOSE2
    for( int ii = 0 ; ii < Ns ; ii ++ )
      printf("%24.16E\n", omegas[ii] );
    #endif
    return omegas;
  }
  

  
}
