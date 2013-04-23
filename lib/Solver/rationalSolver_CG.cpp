/*!
 * @file rationalSolver.cpp
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 * Time-stamp: <2013-04-23 13:10:52 neo>
 */
#include "rationalSolver_CG.hpp"
#include "Communicator/comm_io.hpp"

SolverOutput RationalSolver_CG::solve(Field& sol, const Field& source) const {
  if (Residuals.size() == 0)  {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  Field temp;
  vector_Field shifted_sol;
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

SolverOutput RationalSolver_CG::solve_inv(Field& sol, const Field& source) const {
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  Field temp;
  vector_Field shifted_sol;
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
SolverOutput RationalSolver_CG::solve_noReconstruct(std::vector<Field>& shifted_sol, 
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
SolverOutput RationalSolver_CG::solve_noReconstruct_inv(std::vector<Field>& shifted_sol, 
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

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
SolverOutput RationalSolver_DWF_Optimized::solve(Field& sol, const Field& source) const {
 SolverOutput out;
  if (Residuals.size() == 0)  {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  Field temp;
  vector_Field shifted_sol;
  shifted_sol.resize(Residuals.size());
  sol.resize(source.size());
  
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  temp.resize(source.size());
  
  opr_->solve_ms_eo(shifted_sol, source, out, Poles, Params.MaxIter, Params.GoalPrecision);

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

SolverOutput RationalSolver_DWF_Optimized::solve_inv(Field& sol, const Field& source) const {
  SolverOutput out;
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  Field temp;
  vector_Field shifted_sol;
  shifted_sol.resize(InvResiduals.size());
  sol.resize(source.size());

  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
  temp.resize(source.size());

  opr_->solve_ms_eo(shifted_sol, source, out, InvPoles, Params.MaxIter, Params.GoalPrecision);
 

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
SolverOutput RationalSolver_DWF_Optimized::solve_noReconstruct(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  SolverOutput out;
  if (InvResiduals.size() == 0) {
    CCIO::cout << "[RationalSolver] Inverse Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(InvResiduals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }

  opr_->solve_ms_eo(shifted_sol, source, out, InvPoles, Params.MaxIter, Params.GoalPrecision);

  return out;
}


// Needed in force calculation
// It solves the direct equation (nevertheless the name says "inv" because in the action 
// such a term is associated to the inverse operator
SolverOutput RationalSolver_DWF_Optimized::solve_noReconstruct_inv(std::vector<Field>& shifted_sol, 
						 const Field& source) const {
  SolverOutput out;
  if (Residuals.size() == 0) {
    CCIO::cout << "[RationalSolver] Rational Approximation not initialized yet\n";
    abort();
  }
  shifted_sol.resize(Residuals.size());
 
  for (int i=0; i< shifted_sol.size(); ++i) {
    shifted_sol[i].resize(source.size());
  }
 
  opr_->solve_ms_eo(shifted_sol, source, out, Poles, Params.MaxIter, Params.GoalPrecision);


  return out;
}

#endif
