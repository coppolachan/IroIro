/*!
 * @file solver_CG_DWF_BGQ.cpp
 *
 * @brief Definition of Solver_CG_DWF_Optimized class
 *
 * Time-stamp: <2013-07-04 16:14:45 cossu>
 */

// Valid only for BGQ

#include "solver_CG_DWF_BGQ.hpp"

SolverOutput Solver_CG_DWF_Optimized::solve(Field& xq,const Field& b) const{ 
  SolverOutput Out;
  // direct access to the optimized version inside the Dirac Kernel
  opr_->solve_eo(xq, b, Out, Params.MaxIter,Params.GoalPrecision);
  return Out;
}
