#include "APEsmear.hpp"

using namespace std;

#include "Measurements/GaugeM/staples.h"
#include "include/common_fields.hpp"

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
void Smear_APE::smear(Field& u_smr, const Field& u){

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
//============================================================END=====
