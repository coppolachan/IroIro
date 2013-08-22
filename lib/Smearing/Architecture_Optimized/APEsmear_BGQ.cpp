#include "Smearing/APEsmear.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Geometry/shiftField.hpp"
#include "Geometry/mapping.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"

#include "bgqthread.h"
#include "timings.hpp"
#include "include/messages_macros.hpp"
#include <omp.h>

//====================================================================
std::vector<double> Smear_APE::set_rho(const double common_rho)const{

  std::vector<double> res;

  for(int mn=0; mn<NDIM_*NDIM_; ++mn) res.push_back(common_rho);
  for(int mu=0; mu<NDIM_; ++mu) res[mu + mu*NDIM_] = 0.0;
  return res;
}

//====================================================================
void Smear_APE::smear(GaugeField& u_smr, const GaugeField& u) const{
  using namespace FieldUtils;
  using namespace Mapping;
  double d_rho;
  GaugeField1D Cup;

  u_smr = 0.0;
  init_shiftField();

  GaugeField1D Cdn;
  GaugeField1D WupMu, VupNu;

  double* v_ptr;
  double* w_ptr;
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* Cup_ptr   = Cup.data.getaddr(0);
  double* Cdn_ptr   = Cdn.data.getaddr(0);
  double* u_smr_ptr   = u_smr.data.getaddr(0);
  double* U_ptr = const_cast<GaugeField&>(u).data.getaddr(0);

  int Nvol = CommonPrms::instance()->Nvol();
  const int CC2 = 2*NC_*NC_;

#pragma omp parallel
  { 
    const int nid = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int is = tid*Nvol/nid;
    const int ie = (tid + 1)*Nvol/nid;
    const int ns = ie - is;
    const int str2 = is*CC2;

    for(int mu=0; mu<NDIM_; ++mu){
      v_ptr = U_ptr + mu*Nvol*CC2;
      for(int nu=0; nu<NDIM_; ++nu){
	d_rho = rho[mu + NDIM_ * nu];
	//Explicit staple calculation avoiding temporaries
	w_ptr = U_ptr + nu*Nvol*CC2;
	
	shiftField(WupMu,w_ptr ,mu,Forward());
	shiftField(VupNu,v_ptr ,nu,Forward());
	
	BGWilsonSU3_MatMult_NND(Cup_ptr+str2, w_ptr+str2, VupNu_ptr+str2, WupMu_ptr+str2, ns);
	BGWilsonSU3_MatMult_DNN(VupNu_ptr+str2, w_ptr+str2, v_ptr+str2, WupMu_ptr+str2, ns);
	shiftField(Cdn,VupNu_ptr,nu,Backward());
	BGWilsonSU3_MatAdd(Cup_ptr+str2, 
			   Cdn_ptr+str2,ns); 
	//c++ Cup *= d_rho;
	BGWilsonSU3_MatMultScalar(Cup_ptr+str2, d_rho,ns);

	//c++ AddSlice(u_smr, Cup, mu);
	BGWilsonSU3_MatAdd(u_smr_ptr+str2+Nvol*CC2*mu, 
			   Cup_ptr+str2,ns); 
      }
    }
  }
  
}
//====================================================================
void Smear_APE::derivative(GaugeField& SigmaTerm,
			   const GaugeField& iLambda, 
			   const GaugeField& Gauge) const{
  using namespace Mapping;

  GaugeField1D staple, u_tmp, sh_field ;
  double rho_munu, rho_numu;

  GaugeField1D WupMu, VupNu;
  double* iLambda_mu_ptr;
  double* iLambda_nu_ptr;

  double* U_nu_ptr;
  double* U_mu_ptr;
  double* VupNu_ptr     = VupNu.data.getaddr(0);
  double* WupMu_ptr     = WupMu.data.getaddr(0);
  double* u_tmp_ptr     = u_tmp.data.getaddr(0);
  double* staple_ptr    = staple.data.getaddr(0);
  double* sh_field_ptr  = sh_field.data.getaddr(0);
  double* SigmaTerm_ptr = SigmaTerm.data.getaddr(0);
  double* Gauge_ptr     = const_cast<GaugeField&>(Gauge).data.getaddr(0);
  double* iLambda_ptr   = const_cast<GaugeField&>(iLambda).data.getaddr(0);
 

  int Nvol = CommonPrms::instance()->Nvol();
  const int CC2 = 2*NC_*NC_;

#pragma omp parallel
  { 
    const int nid = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int is = tid*Nvol/nid;
    const int ns = Nvol/nid;
    const int str2 = is*CC2;


    for(int mu = 0; mu < NDIM_; ++mu){
      U_mu_ptr  = Gauge_ptr+ mu*Nvol*CC2;
      iLambda_mu_ptr = iLambda_ptr +mu*Nvol*CC2;
    
      for(int nu = 0; nu < NDIM_; ++nu){
	if(nu==mu) continue;
	U_nu_ptr = Gauge_ptr+ nu*Nvol*CC2;
	iLambda_nu_ptr = iLambda_ptr +nu*Nvol*CC2;
            
	rho_munu = rho[mu + NDIM_ * nu];
	rho_numu = rho[nu + NDIM_ * mu];


	shiftField(WupMu, U_nu_ptr, mu,Forward());
	shiftField(VupNu, U_mu_ptr, nu,Forward());
	BGWilsonSU3_MatMult_NND(staple_ptr+str2, U_nu_ptr+str2, 
				VupNu_ptr+str2, WupMu_ptr+str2, ns);

	BGWilsonSU3_MatMult_DN(u_tmp_ptr+str2, staple_ptr+str2, iLambda_nu_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(u_tmp_ptr+str2, -rho_numu,ns);
	BGWilsonSU3_MatAdd(SigmaTerm_ptr+Nvol*CC2*mu+str2,
			   u_tmp_ptr+str2,ns);

	//-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)

	shiftField(sh_field, iLambda_nu_ptr, mu, Forward());

	BGWilsonSU3_MatMult_ND(u_tmp_ptr+str2, sh_field_ptr+str2, staple_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(u_tmp_ptr+str2, rho_numu,ns);
	BGWilsonSU3_MatAdd(SigmaTerm_ptr+Nvol*CC2*mu+str2,
			   u_tmp_ptr+str2,ns);

	//r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)

	shiftField(sh_field, iLambda_mu_ptr, nu, Forward());

	BGWilsonSU3_MatMult_NND(u_tmp_ptr+str2, U_nu_ptr+str2, sh_field_ptr+str2, U_nu_ptr+str2, ns);
	BGWilsonSU3_MatMult_DN(sh_field_ptr+str2, staple_ptr+str2, u_tmp_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(sh_field_ptr+str2, -rho_munu,ns);
	BGWilsonSU3_MatAdd(SigmaTerm_ptr+Nvol*CC2*mu+str2,
			   sh_field_ptr+str2,ns);

	//-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)

	BGWilsonSU3_MatZero(staple_ptr+str2, ns);
	shiftField(sh_field, U_nu_ptr, mu, Forward());

	BGWilsonSU3_MatMult_NN(u_tmp_ptr+str2, U_mu_ptr+str2,sh_field_ptr+str2, ns);
	BGWilsonSU3_MatMult_DNN(sh_field_ptr+str2, u_tmp_ptr+str2, 
				iLambda_mu_ptr+str2, U_nu_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(sh_field_ptr+str2, -rho_munu,ns);
	BGWilsonSU3_MatAdd(staple_ptr+str2, 
			   sh_field_ptr+str2,ns); 
	BGWilsonSU3_MatMult_DNN(sh_field_ptr+str2, u_tmp_ptr+str2, 
				iLambda_nu_ptr+str2, U_nu_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(sh_field_ptr+str2, rho_numu,ns);
	BGWilsonSU3_MatAdd(staple_ptr+str2, 
			   sh_field_ptr+str2,ns); 

	BGWilsonSU3_MatMult_DN(u_tmp_ptr+str2, U_nu_ptr+str2, iLambda_nu_ptr+str2, ns);

	BGQThread_Barrier(0, nid);

	shiftField(sh_field, u_tmp_ptr, mu, Forward());

	BGWilsonSU3_MatMult_ND(u_tmp_ptr+str2, sh_field_ptr+str2, U_mu_ptr+str2, ns);
	BGWilsonSU3_MatMult_NN(sh_field_ptr+str2, u_tmp_ptr+str2, U_nu_ptr+str2, ns);
	BGWilsonSU3_MatMultScalar(sh_field_ptr+str2,  - rho_numu,ns);
	BGWilsonSU3_MatAdd(staple_ptr+str2, 
			   sh_field_ptr+str2,ns); 


	BGQThread_Barrier(0, nid);
	shiftField(sh_field, staple_ptr, nu, Backward());
	BGWilsonSU3_MatAdd(SigmaTerm_ptr+Nvol*CC2*mu+str2,
			   sh_field_ptr+str2,ns);
      }
    }
  }
}
