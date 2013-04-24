/*!
 * @file solver_Factory.hpp 
 * @brief Declaration of Solver operators factories
 *
 * Time-stamp: <2013-04-23 13:36:36 cossu>
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
class SolverOperatorFactory {
public:
  virtual Solver* getSolver(const Fopr*) = 0;
  virtual Solver* getSolver(const Fopr_Herm*) = 0;//for CG like inverters
  virtual Solver* getSolver(const Fopr_Herm_Precondition*){
    std::cerr<< "getSolver Error: not defined for Fopr_Herm_Precondition Operator" 
	     << std::endl;
    abort();
  }//for Prec like inverters
};

/*!
 * @brief Abstract base class for creating RationalSolver operators
 * It serves as a type definition to distinguish this class of solvers
 */
class RationalSolverOperatorFactory{
public:
  virtual RationalSolver* getSolver(const Fopr*) = 0;
  virtual RationalSolver* getSolver(const Fopr_Herm*) = 0;//for CG like inverters
  virtual RationalSolver* getSolver(const Fopr_Herm_Precondition*){
    std::cerr<< "getSolver Error: not defined for Fopr_Herm_Precondition Operator" 
	     << std::endl;
    abort();
  }//for Prec like inverters
};

/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Conjugate Gradient Solver operator
*/
class SolverCGFactory : public SolverOperatorFactory {
  const XML::node Solver_node;
public:
  SolverCGFactory(const XML::node node):Solver_node(node){}

  Solver_CG* getSolver(const Fopr_Herm* HermitianOperator){
    return new Solver_CG(Solver_node, HermitianOperator);
  }

  Solver_CG* getSolver(const Fopr_Herm_Precondition* HermitianOperator){
    return new Solver_CG(Solver_node, HermitianOperator);
  }

  Solver_CG* getSolver(const Fopr* LinearOperator){
    std::cerr<< "getSolver Error: Solver_CG requires Hermitian Operator" 
	     << std::endl;
    abort();
  }
};

/*!
 @brief Concrete class for creating Conjugate Gradient Solver operator
*/
class SolverCGPrecFactory : public SolverOperatorFactory {
  const XML::node Solver_node;
public:
  SolverCGPrecFactory(const XML::node node):Solver_node(node){}

  Solver_CG_Precondition* getSolver(const Fopr_Herm_Precondition* HermitianOperator){
    return new Solver_CG_Precondition(Solver_node, HermitianOperator);
  }

  Solver_CG_Precondition* getSolver(const Fopr_Herm* HermitianOperator){
    std::cerr<< "getSolver Error: Solver_CG_Precondition requires Hermitian Preconditioned Operator" 
	     << std::endl;
    abort();
  }
  
  Solver_CG_Precondition* getSolver(const Fopr* LinearOperator){
    std::cerr<< "getSolver Error: Solver_CG_Precondition requires Hermitian Preconditioned Operator" 
	     << std::endl;
    abort();
  }
};


#ifdef IBM_BGQ_WILSON
#include "Dirac_ops/dirac_DomainWall_EvenOdd.hpp"
/*!
 @brief Concrete class for creating Conjugate Gradient Solver operator [Optimized version]
*/
class SolverCG_DWF_opt_Factory {
  const XML::node Solver_node;
public:
  SolverCG_DWF_opt_Factory(const XML::node node):Solver_node(node){}

  Solver_CG_DWF_Optimized* getSolver(const Dirac_optimalDomainWall_EvenOdd* DWFopr){
    return new Solver_CG_DWF_Optimized(Solver_node, DWFopr);
  }

};
#endif


/*!
 @brief Concrete class for creating BiConjugate Gradient Stabilized Solver
 operator 
 */
class SolverBiCGStabFactory : public SolverOperatorFactory {
  const XML::node Solver_node;

public:
  SolverBiCGStabFactory(const XML::node node):Solver_node(node){}

  Solver_BiCGStab* getSolver(const Fopr_Herm* HermitianOperator){
    return new Solver_BiCGStab(Solver_node, HermitianOperator);
  }

  Solver_BiCGStab* getSolver(const Fopr_Herm_Precondition* HermitianOperator){
    return new Solver_BiCGStab(Solver_node, HermitianOperator);
  }

  Solver_BiCGStab* getSolver(const Fopr* LinearOperator){
    return new Solver_BiCGStab(Solver_node, LinearOperator);
  }
};

/*!
 @brief Concrete class for creating Rational Conjugate Gradient Solver
 operator 
 */
class RationalSolverCGFactory: public RationalSolverOperatorFactory {
  RaiiFactoryObj<MultiShiftSolver> MS_Solver;

  const XML::node Solver_node;
public:
  RationalSolverCGFactory(const XML::node node):Solver_node(node){}

  RationalSolver_CG* getSolver(const Fopr_Herm* HermitianOperator) {
    MS_Solver.save(new MultiShiftSolver_CG(HermitianOperator,
					   Solver_node));
    return new RationalSolver_CG(MS_Solver.get());
  }

  RationalSolver_CG* getSolver(const Fopr_Herm_Precondition* HermitianOperator) {
    MS_Solver.save(new MultiShiftSolver_CG(HermitianOperator,
					   Solver_node));
    return new RationalSolver_CG(MS_Solver.get());
  }

  RationalSolver_CG*  getSolver(const Fopr* LinearOperator){
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
    return new RationalSolver_DWF_Optimized(DWFopr, Solver_node);
  }

};
#endif

///////////////////////////////////////////////////

namespace SolverOperators{
  SolverOperatorFactory* createSolverOperatorFactory(const XML::node);
  RationalSolverOperatorFactory* createRationalSolverOperatorFactory(const XML::node);
}

#endif 
