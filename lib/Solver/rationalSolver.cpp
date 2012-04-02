/*!
 * @file rationalSolver_CG.cpp
 *
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 */

#include "rationalSolver.hpp"

SolverOutput RationalSolver::solve(Field& sol, const Field& source) const {
  if (Residuals.size() == 0)  {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  SolverOutput out;
  Field temp;
  prop_t shifted_sol;
  shifted_sol.resize(Residuals.size());

  MS_Solver_->solve(shifted_sol, source, Poles, out.diff, out.Iterations);

  // Reconstruct solution (M^dag M)^(a/b)
  sol = source;
  sol *= ConstTerm;

  for (int i = 0; i < Poles.size(); ++i){
    temp = shifted_sol[i];
    temp *= Residuals[i];
    sol += temp;
  }

  return out;
}

SolverOutput RationalSolver::solve_inv(Field& sol, const Field& source) const {
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  SolverOutput out;
  Field temp;
  prop_t shifted_sol;
  shifted_sol.resize(InvResiduals.size());

  MS_Solver_->solve(shifted_sol, source, InvPoles, out.diff, out.Iterations);

  // Reconstruct solution (M^dag M)^(-a/b)
  sol = source;
  sol *= InvConstTerm;

  for (int i = 0; i < InvPoles.size(); ++i){
    temp = shifted_sol[i];
    temp *= InvResiduals[i];
    sol += temp;
  }

  return out;
}

SolverOutput RationalSolver::solve_noReconstruct(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  SolverOutput out;
  Field temp;
  shifted_sol.resize(InvResiduals.size());
  
  MS_Solver_->solve(shifted_sol, source, InvPoles, out.diff, out.Iterations);

  for (int i = 0; i < InvPoles.size(); ++i)
    shifted_sol[i] *= sqrt(InvResiduals[i]);

  return out;
}
