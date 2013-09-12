/*!
 * @file rationalSolver_DWF_BGQ.hpp
 *
 * @brief Wrapper around MultiShiftSolver_CG class
 *
 * Calculates rational functions approximations
 *
 */
#ifndef RATIONALSOLVER_DWF_BGQ_INCLUDED
#define RATIONALSOLVER_DWF_BGQ_INCLUDED

#include "Solver/rationalSolver.hpp"
#include "Solver/multiShiftSolver_CG.hpp"
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"

#ifdef HAVE_LIBBFM
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.hpp"
#endif

class RationalSolver_DWF_Optimized: public RationalSolver {
#ifdef HAVE_LIBBFM
  Dirac_BFM_Wrapper* BFM_opr_;
#endif
  const Dirac_optimalDomainWall_EvenOdd* opr_;
  bool is_BFM;

  int Nsites;
  const MultiShiftSolver_CG_Params Params;

  vector_double Residuals;
  vector_double Poles;
  double ConstTerm;

  vector_double InvResiduals;
  vector_double InvPoles;
  double InvConstTerm;

  void internal_solve(vector_Field&, const Field&, const vector_double&, SolverOutput&) const;

public:
#ifdef HAVE_LIBBFM
  // Standard Constructor
  RationalSolver_DWF_Optimized(Dirac_BFM_Wrapper* BFMopr,
			       const XML::node node);
#endif
  // Standard Constructor
  RationalSolver_DWF_Optimized(const Dirac_optimalDomainWall_EvenOdd* DWFopr,
			       const XML::node node,
			       RationalApprox& RA);

  // To be used in Action constructions
  RationalSolver_DWF_Optimized(const Dirac_optimalDomainWall_EvenOdd* DWFopr,
			       const XML::node node);


  // Destructor
  ~RationalSolver_DWF_Optimized(){};

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
