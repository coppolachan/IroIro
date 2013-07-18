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
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.hpp"

class RationalSolver_DWF_Optimized: public RationalSolver {
  // This duplication of internal pointers must be changed later 
  Dirac_BFM_Wrapper* BFM_opr_;
  bool is_BFM;

  const Dirac_optimalDomainWall_EvenOdd* opr_;
  const MultiShiftSolver_CG_Params Params;

  vector_double Residuals;
  vector_double Poles;
  double ConstTerm;

  vector_double InvResiduals;
  vector_double InvPoles;
  double InvConstTerm;

  void internal_solve(vector_Field&, const Field&, const vector_double&, SolverOutput&) const;

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

  RationalSolver_DWF_Optimized(Dirac_BFM_Wrapper* BFMopr,
			       const XML::node node)
    :BFM_opr_(BFMopr),
     Params(MultiShiftSolver_CG_Params(node)),
     Residuals(0),
     Poles(0),
     ConstTerm(0),
     InvResiduals(0),
     InvPoles(0),
     InvConstTerm(0),
     is_BFM(true){
    BFM_opr_->set_SolverParams(node); 
  };


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
