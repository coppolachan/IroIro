/*!
 * @file dirac_DomainWall_4D.cpp
 *
 * @brief Declaration of Dirac_optimalDomainWall_4D class
 *
 */

#include "dirac_DomainWall_4D.hpp"

const Field Dirac_optimalDomainWall_4D::mult(const Field& f)const{
  int Nconv;
  double diff;
  Field sol5(Dodw_.fsize());
  
  slv_pv_->solve(sol5,Dpv_.mult_dag_prec(Dodw_.mult_prec(Dodw_.Bproj_dag(f))),
                 diff, Nconv);
  return Dodw_.Bproj(sol5);
}

const Field Dirac_optimalDomainWall_4D::mult_dag(const Field& f)const{
  int Nconv;
  double diff;
  Field sol5(Dodw_.fsize());

  slv_pv_->solve(sol5,Dodw_.Bproj_dag(f),diff, Nconv);
  return Dodw_.Bproj(Dodw_.mult_dag_prec(Dpv_.mult(sol5)));
}


const Field Dirac_optimalDomainWall_4D::mult_inv(const Field& f)const{
  int Nconv;
  double diff;
  Field sol5(Dodw_.fsize());
  //prec version to precondition the source
  slv_odw_->solve(sol5,Dodw_.mult_dag_prec(Dpv_.mult_prec(Dodw_.Bproj_dag(f))),
		 diff,Nconv);
  return Dodw_.Bproj(sol5);
}


const Field Dirac_optimalDomainWall_4D::mult_dag_inv(const Field& f)const{
  int Nconv;
  double diff;
  Field sol5(Dodw_.fsize());

  slv_odw_->solve(sol5,Dodw_.Bproj_dag(f),diff,Nconv);
  return Dodw_.Bproj(Dpv_.mult_dag_prec(Dodw_.mult_prec(sol5)));
}

/*!
 * @brief Calculates the sign function of the Kernel
 *
 * It uses the equation
 * \f[{\rm sign}(H_{kernel}) = 2 \gamma_5 D(m=0) - \gamma_5\f]
 *
 */
const Field Dirac_optimalDomainWall_4D::signKernel(const Field& f)const {
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
