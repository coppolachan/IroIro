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

  std::vector<double> Residuals;
  std::vector<double> Poles;
  double ConstTerm;

  std::vector<double> InvResiduals;
  std::vector<double> InvPoles;
  double InvConstTerm;

  RationalSolver(){}; //hide default constructor
public:
  // Standard Constructor
  RationalSolver(const MultiShiftSolver* Solv,
		 RationalApprox& RA)
    :MS_Solver_(Solv),
     Residuals(RA.Residuals()),
     Poles(RA.Poles()),
     ConstTerm(RA.Const()),
     InvResiduals(RA.InvResiduals()),
     InvPoles(RA.InvPoles()),
     InvConstTerm(RA.InvConst()){};

  // Destructor
  ~RationalSolver(){};

  // To be used in Action constructions
  RationalSolver(const MultiShiftSolver* Solv)
    :MS_Solver_(Solv),
     Residuals(0),
     Poles(0),
     ConstTerm(0),
     InvResiduals(0),
     InvPoles(0),
     InvConstTerm(0){};

  void set_Approx(RationalApprox& RA) {
    Residuals = RA.Residuals();
    Poles     = RA.Poles();
    ConstTerm = RA.Const();
    InvResiduals = RA.InvResiduals();
    InvPoles     = RA.InvPoles();
    InvConstTerm = RA.InvConst();
  }

  SolverOutput solve(Field&, const Field&) const;
  SolverOutput solve_inv(Field&, const Field&) const;

  // For force calculations
  SolverOutput solve_noReconstruct(std::vector<Field>&, const Field&) const;
  SolverOutput solve_noReconstruct_inv(std::vector<Field>&, const Field&) const;


  bool check_DdagD() const {return true;}


};




#endif
