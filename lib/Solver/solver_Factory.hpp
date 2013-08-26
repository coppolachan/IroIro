/*!
 * @file solver_Factory.hpp 
 * @brief Declaration of Solver operators factories
 * Time-stamp: <2013-08-21 12:10:08 cossu>
 */
#ifndef SOLVER_FACT_
#define SOLVER_FACT_

#include "Tools/RAIIFactory.hpp"

#include "Solver/solver.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/solver_BiCGStab.hpp"
#include "Solver/multiShiftSolver.hpp"
#include "Solver/rationalSolver.hpp"
#include "Solver/rationalSolver_CG.hpp"
#include "Solver/multiShiftSolver_CG.hpp"

/*!
 * @brief Abstract base class for creating Solver operators
 */
class SolverFactory {
public:
  virtual Solver* getSolver(const Fopr*) = 0;
  virtual Solver* getSolver(const Fopr_Herm*) = 0;//for CG like inverters
};

/*!
 * @brief Abstract base class for creating RationalSolver operators
 * It serves as a type definition to distinguish this class of solvers
 */
class RationalSolverFactory{
public:
  virtual RationalSolver* getSolver(const Fopr*) = 0;
  virtual RationalSolver* getSolver(const Fopr_Herm*) = 0;//for CG like inverters
};

/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Conjugate Gradient Solver operator
*/
class SolverCGFactory : public SolverFactory {
  const XML::node Solver_node;
public:
  SolverCGFactory(const XML::node node):Solver_node(node){}

  Solver_CG* getSolver(const Fopr_Herm* HermitianOp){
    return new Solver_CG(Solver_node, HermitianOp);
  }
   Solver_CG* getSolver(const Fopr* LinearOp){
     std::cerr<< "getSolver Error: Solver_CG requires Hermitian Operator\n";
     abort();
   }
};

#ifdef IBM_BGQ_WILSON
#include "Solver/Architecture_Optimized/solver_CG_DWF_BGQ.hpp"

/*!
 @brief Concrete class for creating Conjugate Gradient Solver operator [Optimized version]
*/
class SolverCG_DWF_opt_Factory {
  const XML::node Solver_node;
public:
  Solver_CG_DWF_Optimized* getSolver(const Dirac_optimalDomainWall_EvenOdd* DWFopr){
    return new Solver_CG_DWF_Optimized(Solver_node,DWFopr);
  }
#ifdef HAVE_LIBBFM
  Solver_CG_DWF_Optimized* getSolver(Dirac_BFM_Wrapper* BFMopr){
    return new Solver_CG_DWF_Optimized(Solver_node,BFMopr);
  }
#endif

  SolverCG_DWF_opt_Factory(const XML::node node):Solver_node(node){};
};
#endif

/*!
 @brief Concrete class for creating BiConjugate Gradient Stabilized Solver
 operator 
 */
class SolverBiCGStabFactory : public SolverFactory {
  const XML::node Solver_node;

public:
  SolverBiCGStabFactory(const XML::node node):Solver_node(node){}

  Solver_BiCGStab* getSolver(const Fopr_Herm* HermitianOp){
    return new Solver_BiCGStab(Solver_node,HermitianOp);
  }
  Solver_BiCGStab* getSolver(const Fopr* LinearOp){
    return new Solver_BiCGStab(Solver_node, LinearOp);
  }
};

/*!
 @brief Concrete class for creating Rational Conjugate Gradient Solver
 operator 
 */
class RationalSolverCGFactory: public RationalSolverFactory {
  RaiiFactoryObj<MultiShiftSolver> MS_Solver;

  const XML::node Solver_node;
public:
  RationalSolverCGFactory(const XML::node node):Solver_node(node){}

  RationalSolver_CG* getSolver(const Fopr_Herm* HermitianOp) {
    MS_Solver.save(new MultiShiftSolver_CG(HermitianOp,Solver_node));
    return new RationalSolver_CG(MS_Solver.get());
  }
  RationalSolver_CG*  getSolver(const Fopr* LinearOp){
    std::cerr<< "getSolver Error: RationalSolver requires Hermitian Operator" 
	     << std::endl;
    abort();
  }
};

/*!
 @brief Concrete class for creating Rational Conjugate Gradient Solver
 operator optimized version for DWF on BGQ
 */
#ifdef IBM_BGQ_WILSON
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
#include "Solver/Architecture_Optimized/rationalSolver_DWF_BGQ.hpp"

class RationalSolverCGFactory_DWF_Optimized {
   const XML::node Solver_node;
public:
  RationalSolverCGFactory_DWF_Optimized(const XML::node node):Solver_node(node){}

  RationalSolver_DWF_Optimized* getSolver(const Dirac_optimalDomainWall_EvenOdd* DWFopr) {
    return new RationalSolver_DWF_Optimized(DWFopr,Solver_node);
  }
#ifdef HAVE_LIBBFM
  RationalSolver_DWF_Optimized* getSolver(Dirac_BFM_Wrapper* BFMopr){
    return new RationalSolver_DWF_Optimized(BFMopr, Solver_node);
  }
#endif

};
#endif

///////////////////////////////////////////////////

namespace Solvers{
  SolverFactory* createSolverFactory(const XML::node);
  RationalSolverFactory* createRationalSolverFactory(const XML::node);
}

#endif 
