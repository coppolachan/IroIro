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
  SolverOutput out;
  prop_t shifted_sol;
  shifted_sol.resize(Residuals.size());
  for (int i = 0; i < Residuals.size(); ++i)
    shifted_sol.push_back(sol);

  MS_Solver_->solve(shifted_sol, source, Poles, out.diff, out.Iterations);

  //shifted_sol contains the solutions of (M - Poles[i])x = source

  // Reconstruct solution (M^dag M)^(a/b)
  Field temp;

  sol = source;
  sol *= ConstTerm;

  for (int i = 0; i < Poles.size(); ++i){
    temp = shifted_sol[i];
    temp *= Residuals[i];
    sol += temp;
  }

  return out;
}
