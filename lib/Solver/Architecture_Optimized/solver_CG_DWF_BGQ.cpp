/*!
 * @file solver_CG_DWF_BGQ.cpp
 *
 * @brief Definition of Solver_CG_DWF_Optimized class
 *
 * Time-stamp: <2013-07-16 10:57:15 cossu>
 */

// Valid only for BGQ

#include "solver_CG_DWF_BGQ.hpp"

SolverOutput Solver_CG_DWF_Optimized::solve(Field& xq,const Field& b) const{ 
  SolverOutput Out;
  // direct access to the optimized version inside the Dirac Kernel
  if (is_BFM) {
    FermionField xq_FF(xq);
    FermionField b_FF(b);
    BFM_opr_->solve_CGNE(xq_FF, b_FF);
    xq = xq_FF.data;
  } else {
  opr_->solve_eo(xq, b, Out, Params.MaxIter,Params.GoalPrecision);
  }
  return Out;
}
