#include "APEsmear.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "Main/Geometry/mapper.hpp"
#include "Tools/sunMatUtils.hpp"


//====================================================================
std::vector<double> Smear_APE::set_rho(const double common_rho)const{

  std::vector<double> res(NDIM_*NDIM_);

  for(int mu = 0; mu < NDIM_; ++mu){
   for(int nu = 0; nu < NDIM_; ++nu){
     res[mu + nu*NDIM_] = common_rho;
   }
  }
  for(int mu = 0; mu < NDIM_; ++mu){
     res[mu + mu*NDIM_] = 0.0;
  }

  return res;
}

//====================================================================
void Smear_APE::smear(GaugeField& u_smr, const GaugeField& u) const{
  using namespace FieldUtils;
  double d_rho;
  GaugeField1D Cup, Cdn;
  Staples stpl;

  u_smr = 0.0;

  for(int mu = 0; mu < NDIM_; ++mu){
   for(int nu = 0; nu < NDIM_; ++nu){
     Cup = stpl.upper(u,mu,nu);
     Cdn = stpl.lower(u,mu,nu);

     d_rho = rho[mu + NDIM_ * nu];

     Cup += Cdn;
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
  using namespace MapsEnv;

  Staples stpl;
  GaugeField1D staple, u_tmp, iLambda_mu, iLambda_nu;
  GaugeField1D U_mu, U_nu;
  GaugeField1D sh_field;
  SUNmat temp_mat, temp_mat2;
  double rho_munu, rho_numu;

  int Nvol = CommonPrms::instance()->Nvol();

  for(int mu = 0; mu < NDIM_; ++mu){
    //FieldUP UpMu(&staple.Format,mu);//1D field shifter
    U_mu       = DirSlice(  Gauge, mu);
    iLambda_mu = DirSlice(iLambda, mu);
    
    for(int nu = 0; nu < NDIM_; ++nu){
      if(nu==mu) continue;
      
      //FieldDN DnNu(&staple.Format,nu);
      //FieldUP UpNu(&staple.Format,nu);
      U_nu       = DirSlice(  Gauge, nu);
      iLambda_nu = DirSlice(iLambda, nu);
      
      rho_munu = rho[mu + NDIM_ * nu];
      rho_numu = rho[nu + NDIM_ * mu];
      
      staple = stpl.upper(Gauge,mu,nu);
      
      for (int site = 0; site < Nvol; ++site){
	temp_mat = matrix_dag(staple,site) * matrix(iLambda_nu,site);
	temp_mat *= - rho_numu;
	AddMatrix(SigmaTerm, temp_mat, site, mu);
      }//-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)

      sh_field = shift(iLambda_nu, mu, Forward);
      for (int site = 0; site < Nvol; ++site){
	temp_mat = matrix(sh_field,site) * matrix_dag(staple,site);
	temp_mat *= rho_numu;
	AddMatrix(SigmaTerm, temp_mat, site, mu);
      }//r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
 
      sh_field = shift(iLambda_mu, nu, Forward);
      for (int site = 0; site < Nvol; ++site){
	temp_mat = matrix(U_nu,site) * matrix(sh_field,site) * matrix_dag(U_nu,site);
	temp_mat = matrix_dag(staple,site) * temp_mat;
	temp_mat *= - rho_munu;
	AddMatrix(SigmaTerm, temp_mat, site, mu);
      }//-r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)
      
      
      staple = 0.0;
      sh_field = shift(U_nu, mu, Forward);
      for (int site = 0; site < Nvol; ++site){
	temp_mat2 = matrix_dag(sh_field,site) * matrix_dag(U_mu,site);
	temp_mat = temp_mat2  * matrix(iLambda_mu,site) * matrix(U_nu,site);
	temp_mat *= - rho_munu;
	AddMatrix(staple, temp_mat, site);
	temp_mat = temp_mat2 * matrix(iLambda_nu,site) * matrix(U_nu,site);
	temp_mat *= rho_numu;
	AddMatrix(staple, temp_mat, site);
      }     

      for (int site = 0; site < Nvol; ++site){
	temp_mat = matrix_dag(U_nu,site) * matrix(iLambda_nu,site);
	SetMatrix(u_tmp, temp_mat, site);
      }     
      
      sh_field = shift(u_tmp, mu, Forward);
      for (int site = 0; site < Nvol; ++site){
	temp_mat = matrix(sh_field,site) * matrix_dag(U_mu,site) * matrix(U_nu,site);
	temp_mat *= - rho_numu;
	AddMatrix(staple, temp_mat, site);
      }     

      sh_field = shift(staple, nu, Backward);
      
      AddSlice(SigmaTerm, sh_field, mu);

    }
  }
}



//====================================================================
//============================================================END=====
