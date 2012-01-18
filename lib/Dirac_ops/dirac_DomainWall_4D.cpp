/*!
 * @file dirac_DomainWall_4D.cpp
 *
 * @brief Declaration of Dirac_optimalDomainWall_4D class
 *
 */

#include "dirac_DomainWall_4D.hpp"

const Field Dirac_optimalDomainWall_4D::mult(const Field& f)const{
  // D_dw(m)
  Field t5(Dodw_.fsize());
  t5 = Dodw_.mult(Dodw_.Bproj_dag(f));
  // Dpv_^-1
  Field src(Dodw_.fsize());
  src = Dpv_.mult_dag_prec(Dpv_.left_prec(t5));
  Field sol5(Dodw_.fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
  t5 = Dpv_.right_prec(sol5);
#if VERBOSITY > 0
  monitor.print();
#endif
  return Dodw_.Bproj(t5);
}

const Field Dirac_optimalDomainWall_4D::mult_dag(const Field& f)const{
  // Dpv_^dag^-1
  Field src(Dodw_.fsize());
  src = Dpv_.right_dag_prec(Dodw_.Bproj_dag(f));
  Field sol5(Dodw_.fsize());
  SolverOutput monitor = slv_pv_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  Field t5(Dodw_.fsize());
  t5 = Dpv_.left_dag_prec(Dpv_.mult_prec(sol5));
  // D_dw^dag(m)
  return Dodw_.Bproj(Dodw_.mult_dag(t5));
}

const Field Dirac_optimalDomainWall_4D::mult_inv(const Field& f)const{
  Field t5(Dodw_.fsize());
  // D_pv
  t5 = Dpv_.mult(Dodw_.Bproj_dag(f));
  // D_dw^-1
  Field src(Dodw_.fsize());
  src = Dodw_.mult_dag_prec(Dodw_.left_prec(t5));
  Field sol5(Dodw_.fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  t5 = Dodw_.right_prec(sol5);
  return Dodw_.Bproj(t5);
}

const Field Dirac_optimalDomainWall_4D::mult_dag_inv(const Field& f)const{
  // D_dw^dag^-1
  Field src(Dodw_.fsize());
  src = Dodw_.right_dag_prec(Dodw_.Bproj_dag(f));
  Field sol5(Dodw_.fsize());
  SolverOutput monitor = slv_odw_->solve(sol5,src);
#if VERBOSITY > 0
  monitor.print();
#endif
  Field t5(Dodw_.fsize());
  t5 = Dodw_.left_dag_prec(Dodw_.mult_prec(sol5));
  // D_pv^dag
  return Dodw_.Bproj(Dpv_.mult_dag(t5));
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
