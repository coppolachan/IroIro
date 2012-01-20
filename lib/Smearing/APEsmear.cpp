#include "APEsmear.hpp"

using namespace std;

#include "Measurements/GaugeM/staples.h"
#include "include/common_fields.hpp"

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

     u_smr += Cup.U;

   }
  }
}
//====================================================================
//============================================================END=====
