/*!
 * @file rationalSolver_CG.hpp
 *
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 */

#ifndef RATIONALSOLVER_CG_INCLUDED
#define RATIONALSOLVER_CG_INCLUDED

#include "rationalSolver.hpp"

class RationalSolver_CG: public RationalSolver {
  const MultiShiftSolver* MS_Solver_;

  vector_double Residuals;
  vector_double Poles;
  double ConstTerm;

  vector_double InvResiduals;
  vector_double InvPoles;
  double InvConstTerm;

  RationalSolver_CG(){}; //hide default constructor
public:
  // Standard Constructor
  RationalSolver_CG(const MultiShiftSolver* Solv,
		 RationalApprox& RA)
    :MS_Solver_(Solv),
     Residuals(RA.Residuals()),
     Poles(RA.Poles()),
     ConstTerm(RA.Const()),
     InvResiduals(RA.InvResiduals()),
     InvPoles(RA.InvPoles()),
     InvConstTerm(RA.InvConst()){};

  // Destructor
  ~RationalSolver_CG(){};

  // To be used in Action constructions
  RationalSolver_CG(const MultiShiftSolver* Solv)
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
  SolverOutput solve_noReconstruct(vector_Field&, const Field&) const;
  SolverOutput solve_noReconstruct_inv(vector_Field&, const Field&) const;

  bool check_DdagD() const {return true;}
};

#ifdef IBM_BGQ_WILSON
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
class RationalSolver_DWF_Optimized: public RationalSolver {
  const Dirac_optimalDomainWall_EvenOdd* opr_;
  const MultiShiftSolver_CG_Params Params;

  vector_double Residuals;
  vector_double Poles;
  double ConstTerm;

  vector_double InvResiduals;
  vector_double InvPoles;
  double InvConstTerm;

public:
  // Standard Constructor
  RationalSolver_DWF_Optimized(const Dirac_optimalDomainWall_EvenOdd* DWFopr,
			       const XML::node node,
			       RationalApprox& RA)
    :opr_(DWFopr),
     Params(MultiShiftSolver_CG_Params(node)),
     Residuals(RA.Residuals()),
     Poles(RA.Poles()),
     ConstTerm(RA.Const()),
     InvResiduals(RA.InvResiduals()),
     InvPoles(RA.InvPoles()),
     InvConstTerm(RA.InvConst()){};

  // Destructor
  ~RationalSolver_DWF_Optimized(){};

  // To be used in Action constructions
  RationalSolver_DWF_Optimized(const Dirac_optimalDomainWall_EvenOdd* DWFopr,
			       const XML::node node)
    :opr_(DWFopr),
     Params(MultiShiftSolver_CG_Params(node)),
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
  SolverOutput solve_noReconstruct(vector_Field&, const Field&) const;
  SolverOutput solve_noReconstruct_inv(vector_Field&, const Field&) const;


  bool check_DdagD() const {return true;}


};

#endif

#endif

