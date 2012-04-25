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
  Field temp;
  prop_t shifted_sol;
  shifted_sol.resize(Residuals.size());
  sol.resize(source.size());
  
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  temp.resize(source.size());
  
  SolverOutput out = MS_Solver_->solve(shifted_sol, source, Poles, out.diff, out.Iterations);

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
  Field temp;
  prop_t shifted_sol;
  shifted_sol.resize(InvResiduals.size());
  sol.resize(source.size());

  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  temp.resize(source.size());


  SolverOutput out = MS_Solver_->solve(shifted_sol, source, InvPoles, out.diff, out.Iterations);

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


// Needed in force calculation
// It solves the inverse equation
SolverOutput RationalSolver::solve_noReconstruct(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(InvResiduals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
 
  SolverOutput out = MS_Solver_->solve(shifted_sol, source, InvPoles, out.diff, out.Iterations);

  return out;
}


// Needed in force calculation
// It solves the direct equation (nevertheless the name says "inv" because in the action 
// such a term is associated to the inverse operator
SolverOutput RationalSolver::solve_noReconstruct_inv(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  if (Residuals.size() == 0) {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(Residuals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
 
  SolverOutput out = MS_Solver_->solve(shifted_sol, source, Poles, out.diff, out.Iterations);

  return out;
}

