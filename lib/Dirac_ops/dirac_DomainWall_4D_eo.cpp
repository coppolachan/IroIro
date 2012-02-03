/*!
 * @file dirac_DomainWall_4D_eo.cpp
 * @brief Methods of Dirac_optimalDomainWall_4D_eo class
 */
#include "dirac_DomainWall_4D_eo.hpp"

const Field Dirac_optimalDomainWall_4D_eo::mult(const Field& f)const{
  // Dpv_^-1
  Field sol5(f5size_);
  invDpv_.invert(sol5,D_.mult(D_.Bproj_dag(f)));
  return D_.Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eo::mult_dag(const Field& f)const{
  // Dpv_^dag^-1
  Field sol5(f5size_);
  invDpv_.invert_dag(sol5,D_.mult_dag(D_.Bproj_dag(f)));
  return D_.Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eo::mult_inv(const Field& f)const{
  // Dodw_^-1
  Field sol5(f5size_);
  invD_.invert(sol5,Dpv_.mult(D_.Bproj_dag(f)));
  return D_.Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D_eo::mult_dag_inv(const Field& f)const{
  // D_odw^dag^-1
  Field sol5(f5size_);
  invD_.invert_dag(sol5,Dpv_.mult(D_.Bproj_dag(f)));
  return D_.Bproj(sol5);
}

/*!
 * @brief Calculates the sign function of the Kernel
 *
 * It uses the equation
 * \f[{\rm sign}(H_{kernel}) = 2 \gamma_5 D(m=0) - \gamma_5\f]
 */
const Field Dirac_optimalDomainWall_4D_eo::signKernel(const Field& f)const {
  double mass = Dodw_.getMass();
  Field signK = mult(f);
  Field aux = f;
  aux *= mass;

  signK -= aux;  //signK = D(m) - m 
  signK /= (1.0 - mass);//signK = (D(m) - m)/(1-m) = D(0)

  signK = gamma5(signK); 
  signK *= 2.0; // 2 \gamma_5 D(m=0) 

  aux = gamma5(f);
  signK -= aux;
  
  return signK;
}
