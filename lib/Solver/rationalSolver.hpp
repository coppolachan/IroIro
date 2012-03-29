/*!
 * @file rationalSolver_CG.hpp
 *
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 */

#ifndef RATIONALSOLVER_INCLUDED
#define RATIONALSOLVER_INCLUDED

#include <vector>
#include "solver.hpp"
#include "multiShiftSolver_CG.hpp"
#include "Tools/RationalApprox/rationalapprox.hpp"

class RationalSolver: public Solver {
  const MultiShiftSolver* MS_Solver_;
  const RationalApprox& RApprox_;

  const std::vector<double> Residuals;
  const std::vector<double> Poles;
  const double ConstTerm;

public:
  // Constructor
  RationalSolver(const MultiShiftSolver* Solv,
		 const RationalApprox& RA)
    :MS_Solver_(Solv),
     RApprox_(RA),
     Residuals(RApprox_.Residuals()),
     Poles(RApprox_.Poles()),
     ConstTerm(RApprox_.Const()){};
  ~RationalSolver(){};

  SolverOutput solve(Field&, const Field&) const;
  bool check_DdagD() const {return true;}


};




#endif
