/*!
 * @file solver_CG_DWF_BGQ.hpp
 *
 * @brief Declaration of Solver_CG_DWF_Optimized class
 *
 * These Solvers have direct access to the Dirac Operator class
 * in order to use special optimized functions without the
 * level of the Fopr
 *
 * Time-stamp: <2013-07-04 16:10:20 cossu>
 */
#ifndef SOLVER_CG_BGQ_INCLUDED
#define SOLVER_CG_BGQ_INCLUDED

#include "include/pugi_interface.h"
#include "Solver/solver.hpp"
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
#include "Dirac_ops/BFM_Wrapper/dirac_BFM_wrapper.hpp"
#include "Solver/solver_CG_params.hpp"
/*!
 * @brief Solves \f$Dx = b\f$ using 
 * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient_method">Conjugate Gradient method</a>
 *
 * Domain Wall fermions optimized operator for BGQ
 * Wrapper class
 * 
 * WARNING: Never use with (non trivially) preconditioned operators
 */
class Solver_CG_DWF_Optimized: public Solver{
private:
  // This duplication of internal pointers must be changed later 
  Dirac_BFM_Wrapper* BFM_opr_;
  bool is_BFM;

  const Dirac_optimalDomainWall_EvenOdd* opr_;
  const Solver_CG_Prms Params;/*!< @brief Inputs container */

public:
  // Create different constructors for BFM
  // pass a BFM dirac operator
  Solver_CG_DWF_Optimized(const XML::node Solver_node,
			  Dirac_BFM_Wrapper* BFMKernel)
    :BFM_opr_(BFMKernel),
     Params(Solver_CG_Prms(Solver_node)),
     is_BFM(true){
    BFM_opr_->set_SolverParams(Solver_node);
      }
  

  Solver_CG_DWF_Optimized(const double prec,const int MaxIterations,
			  const Dirac_optimalDomainWall_EvenOdd* DWFopr)
    :opr_(DWFopr),
     Params(Solver_CG_Prms(prec, MaxIterations)),
     is_BFM(false){}

  Solver_CG_DWF_Optimized(const XML::node Solver_node,
			  const Dirac_optimalDomainWall_EvenOdd* DWFopr)
    :opr_(DWFopr),
     Params(Solver_CG_Prms(Solver_node)),
     is_BFM(false){}

  ~Solver_CG_DWF_Optimized(){}

  SolverOutput solve(Field& solution, const Field& source) const;

  bool check_DdagD() const {
    //empty in this case
    return 1;
  }
};

#endif



