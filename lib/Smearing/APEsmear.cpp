#include "APEsmear.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"

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
  double d_rho;
  GaugeField1D Cup, Cdn;
  Staples stpl;

  u_smr = 0.0;
  for(int mu=0; mu<NDIM_; ++mu){
    for(int nu=0; nu<NDIM_; ++nu){

      //Cup = stpl.upper(u,mu,nu);
      //Cdn = stpl.lower(u,mu,nu);
      //Cup += Cdn;
      Cup = stpl.upper_lower(u,mu,nu);

      d_rho = rho[mu + NDIM_ * nu];
      Cup *= d_rho;
      AddSlice(u_smr, Cup, mu);
    }
  }
}
//====================================================================
void Smear_APE::derivative(GaugeField& SigmaTerm,
			   const GaugeField& iLambda, 
			   const GaugeField& Gauge) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  using namespace Mapping;

  Staples stpl;
  GaugeField1D staple, u_tmp, iLambda_mu, iLambda_nu;
  GaugeField1D U_mu, U_nu;
  GaugeField1D sh_field;
  SUNmat temp_mat, temp_mat2;
  double rho_munu, rho_numu;

#ifdef IBM_BGQ_WILSON
  double* iLambda_mu_ptr = iLambda_mu.data.getaddr(0);
  double* iLambda_nu_ptr = iLambda_nu.data.getaddr(0);
  double* u_tmp_ptr      = u_tmp.data.getaddr(0);
  double* staple_ptr     = staple.data.getaddr(0);
  double* sh_field_ptr   = sh_field.data.getaddr(0);
  double* U_nu_ptr       = U_nu.data.getaddr(0);
  double* U_mu_ptr       = U_mu.data.getaddr(0);
#endif

  int Nvol = CommonPrms::instance()->Nvol();

  for(int mu = 0; mu < NDIM_; ++mu){
    U_mu       = DirSlice(  Gauge, mu);
    iLambda_mu = DirSlice(iLambda, mu);
    
    for(int nu = 0; nu < NDIM_; ++nu){
      if(nu==mu) continue;
      
      U_nu       = DirSlice(  Gauge, nu);
      iLambda_nu = DirSlice(iLambda, nu);
      
      rho_munu = rho[mu + NDIM_ * nu];
      rho_numu = rho[nu + NDIM_ * mu];
      
      staple = stpl.upper(Gauge,mu,nu);
      
#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_DN(u_tmp_ptr, staple_ptr, iLambda_nu_ptr, Nvol);
      u_tmp *= -rho_numu;
      AddSlice(SigmaTerm, u_tmp, mu);
#else
      for (int site = 0; site < Nvol; ++site){
	temp_mat = mat_dag(staple,site) * mat(iLambda_nu,site);
	temp_mat *= - rho_numu;
	AddMat(SigmaTerm, temp_mat, site, mu);
      }//-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
#endif

      sh_field = shiftField(iLambda_nu, mu, Forward());

#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_ND(u_tmp_ptr, sh_field_ptr, staple_ptr, Nvol);
      u_tmp *= rho_numu;
      AddSlice(SigmaTerm, u_tmp, mu);
#else
      for (int site = 0; site < Nvol; ++site){
	temp_mat = mat(sh_field,site) * mat_dag(staple,site);
	temp_mat *= rho_numu;
	AddMat(SigmaTerm, temp_mat, site, mu);
      }//r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
#endif
 
      sh_field = shiftField(iLambda_mu, nu, Forward());

#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_NND(u_tmp_ptr, U_nu_ptr, sh_field_ptr, U_nu_ptr, Nvol);
      BGWilsonSU3_MatMult_DN(sh_field_ptr, staple_ptr, u_tmp_ptr, Nvol);
      sh_field *= -rho_munu;
      AddSlice(SigmaTerm, sh_field, mu);
#else
      for (int site = 0; site < Nvol; ++site){
	temp_mat = mat(U_nu,site) * mat(sh_field,site) * mat_dag(U_nu,site);
	temp_mat = mat_dag(staple,site) * temp_mat;
	temp_mat *= - rho_munu;
	AddMat(SigmaTerm, temp_mat, site, mu);
      }//-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
#endif
      
      staple = 0.0;
      sh_field = shiftField(U_nu, mu, Forward());

#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_NN(u_tmp_ptr, U_mu_ptr,sh_field_ptr, Nvol);
      BGWilsonSU3_MatMult_DNN(sh_field_ptr, u_tmp_ptr, iLambda_mu_ptr, U_nu_ptr, Nvol);
      sh_field *= -rho_munu;
      staple += sh_field;
      BGWilsonSU3_MatMult_DNN(sh_field_ptr, u_tmp_ptr, iLambda_nu_ptr, U_nu_ptr, Nvol);
      sh_field *= rho_numu;
      staple += sh_field;
#else
      for (int site = 0; site < Nvol; ++site){
	temp_mat2 = mat_dag(sh_field,site) * mat_dag(U_mu,site);
	temp_mat = temp_mat2  * mat(iLambda_mu,site) * mat(U_nu,site);
	temp_mat *= - rho_munu;
	AddMat(staple, temp_mat, site);
	temp_mat = temp_mat2 * mat(iLambda_nu,site) * mat(U_nu,site);
	temp_mat *= rho_numu;
	AddMat(staple, temp_mat, site);
      } 
#endif 

#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_DN(u_tmp_ptr, U_nu_ptr, iLambda_nu_ptr, Nvol);
#else   
      for (int site = 0; site < Nvol; ++site){
	temp_mat = mat_dag(U_nu,site) * mat(iLambda_nu,site);
	SetMat(u_tmp, temp_mat, site);
      } 
#endif
   
      sh_field = shiftField(u_tmp, mu, Forward());
#ifdef IBM_BGQ_WILSON
      BGWilsonSU3_MatMult_ND(u_tmp_ptr, sh_field_ptr, U_mu_ptr, Nvol);
      BGWilsonSU3_MatMult_NN(sh_field_ptr, u_tmp_ptr, U_nu_ptr, Nvol);
      sh_field *= - rho_numu;
      staple += sh_field;
#else   
      for (int site = 0; site < Nvol; ++site){
	temp_mat = mat(sh_field,site) * mat_dag(U_mu,site) * mat(U_nu,site);
	temp_mat *= - rho_numu;
	AddMat(staple, temp_mat, site);
      }   
#endif
  
      sh_field = shiftField(staple, nu, Backward());
      AddSlice(SigmaTerm, sh_field, mu);
    }
  }
}
