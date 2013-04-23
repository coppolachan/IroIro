/*!
 * @file rationalSolver.hpp
 *
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 */

#ifndef RATIONALSOLVER_INCLUDED
#define RATIONALSOLVER_INCLUDED

#include "solver.hpp"
#include "multiShiftSolver_CG.hpp"
#include "Tools/RationalApprox/rationalapprox.hpp"

class RationalSolver: public Solver {
public:
  virtual SolverOutput solve_inv(Field&,
				 const Field&) const = 0;

  // For force calculations
  virtual SolverOutput solve_noReconstruct(vector_Field&,
					   const Field&) const = 0;
  virtual SolverOutput solve_noReconstruct_inv(vector_Field&, 
					       const Field&) const = 0;
  virtual void set_Approx(RationalApprox& RA) = 0;
};



#endif
