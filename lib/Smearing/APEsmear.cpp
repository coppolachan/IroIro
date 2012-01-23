#include "APEsmear.hpp"

using namespace std;

#include "Measurements/GaugeM/staples.h"
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

typedef ShiftField_up<GaugeFieldFormat> FieldUP;
typedef ShiftField_dn<GaugeFieldFormat> FieldDN;

//====================================================================
std::valarray<double> Smear_APE::set_rho(const double common_rho)const{

  std::valarray<double> res(Ndim*Ndim);

  for(int mu = 0; mu < Ndim; ++mu){
   for(int nu = 0; nu < Ndim; ++nu){
     res[mu + nu*Ndim] = common_rho;
   }
  }
  for(int mu = 0; mu < Ndim; ++mu){
     res[mu + mu*Ndim] = 0.0;
  }

  return res;
}

//====================================================================
void Smear_APE::smear(Field& u_smr, const Field& u) const{

  double d_rho;
  GaugeField1D Cup, Cdn;
  Staples stpl(Gformat);

  u_smr = 0.0;

  for(int mu = 0; mu < Ndim; ++mu){
   for(int nu = 0; nu < Ndim; ++nu){

     Cup.U = stpl.upper(u,mu,nu);
     Cdn.U = stpl.lower(u,mu,nu);

     d_rho = rho[mu + Ndim * nu];


     Cup.U += Cdn.U;
     Cup.U *= d_rho;
     u_smr.add(Gformat.dir_slice(mu),Cup.U.getva());
 
   }
  }
}
//====================================================================
void Smear_APE::derivative(Field& SigmaTerm,
			   const Field& iLambda, 
			   const Field& Gauge) const{
  using namespace SUNmat_utils;

  Staples stpl(Gformat);
  GaugeField1D staple, u_tmp, iLambda_mu, iLambda_nu;
  GaugeField1D U_mu, U_nu;
  SUNmat temp_mat, temp_mat2;
  double rho_munu, rho_numu;

  int Nvol = CommonPrms::instance()->Nvol();

  for(int mu = 0; mu < Ndim; ++mu){
    FieldUP UpMu(&staple.Format,mu);
    U_mu = Gauge[Gformat.dir_slice(mu)]; 
    iLambda_mu = iLambda[Gformat.dir_slice(mu)];
    
    for(int nu = 0; nu < Ndim; ++nu){
      if(nu==mu) continue;
      
      FieldDN DnNu(&staple.Format,nu);
      FieldUP UpNu(&staple.Format,nu);
      U_nu = Gauge[Gformat.dir_slice(nu)];
      iLambda_nu = iLambda[Gformat.dir_slice(nu)];
      
      rho_munu = rho[mu + Ndim * nu];
      rho_numu = rho[nu + Ndim * mu];
      
      staple.U = stpl.upper(Gauge,mu,nu);
      
      for (int site = 0; site < Nvol; ++site){
	temp_mat = u_dag(staple,site) * u(iLambda_nu,site);
	temp_mat *= - rho_numu;
	SigmaTerm.add(Gformat.cslice(0,site,mu),temp_mat.getva());
      }//-r_numu*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)*Lambda_nu(x)
      
      
      UpMu.setf(iLambda_nu.U);
      for (int site = 0; site < Nvol; ++site){
	temp_mat = u(UpMu,site) * u_dag(staple,site);
	temp_mat *= rho_numu;
	SigmaTerm.add(Gformat.cslice(0,site,mu),temp_mat.getva());
      }//r_numu*Lambda_nu(mu)*U_nu(x+mu)*Udag_mu(x+nu)*Udag_nu(x)
      
      
      UpNu.setf(iLambda_mu.U);
      for (int site = 0; site < Nvol; ++site){
	temp_mat = u(U_nu,site) * u(UpNu,site) * u_dag(U_nu,site);
	temp_mat = u_dag(staple,site) * temp_mat;
	temp_mat *= - rho_munu;
	SigmaTerm.add(Gformat.cslice(0,site,mu),temp_mat.getva());
      }//r_munu*U_nu(x+mu)*Udag_mu(x+nu)*Lambda_mu(x+nu)*Udag_nu(x)

      staple.U = 0.0;
      UpMu.setf(U_nu.U);
      for (int site = 0; site < Nvol; ++site){
	temp_mat2 = u_dag(UpNu,site) * u_dag(U_mu,site);
	temp_mat = temp_mat2  * u(iLambda_mu,site) * u(U_nu,site);
	temp_mat *= - rho_munu;
	staple.U.add(staple.Format.cslice(0,site,0),temp_mat.getva());
	temp_mat = temp_mat2 * u(iLambda_nu,site) * u(U_nu,site);
	temp_mat *= rho_numu;
	staple.U.add(staple.Format.cslice(0,site,0),temp_mat.getva());
      }     

      for (int site = 0; site < Nvol; ++site){
	temp_mat = u_dag(U_nu,site) * u(iLambda_nu,site);
	u_tmp.U.set(u_tmp.Format.cslice(0,site,0),temp_mat.getva());
      }     
      UpMu.setf(u_tmp.U);

      for (int site = 0; site < Nvol; ++site){
	temp_mat = u(u_tmp,site) * u_dag(U_mu,site) * u(U_nu,site);
	temp_mat *= - rho_numu;
	staple.U.add(staple.Format.cslice(0,site,0),temp_mat.getva());
      }     

      DnNu.setf(staple.U);
      
      SigmaTerm.add(Gformat.dir_slice(mu),DnNu.getva());

    }
  }
}



//====================================================================
//============================================================END=====
