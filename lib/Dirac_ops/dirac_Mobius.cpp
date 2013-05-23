/*! @file dirac_Mobius.cpp
 * @brief implementation of the Dirac_Mobius class 
 Time-stamp: <2013-05-21 11:52:03 noaki> 
 */
#include "dirac_Mobius.hpp"
#include "Fields/field_expressions.hpp"
#include "Solver/solver.hpp"

const Field Dirac_Mobius::mult(const Field& f) const{
  using namespace FieldExpression;
  Field src = Dd_->mult_dag(f);
  Field w(Dd_->fsize());
  SolverOutput monitor = slv_->solve(w,src);
  w *= c2_;
  w += c1_*f;
  return w;
}

const Field Dirac_Mobius::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

