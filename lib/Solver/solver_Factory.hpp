/*!
 * @file solver_Factory.hpp 
 *
 * @brief Declaration of Solver operators factories
 */

#ifndef SOLVER_FACT_
#define SOLVER_FACT_

#include "Tools/RAIIFactory.hpp"

#include "Solver/solver.h"
#include "Solver/solver_CG.h"
#include "Solver/solver_BiCGStab.h"

/*!
 * @brief Abstract base class for creating Solver operators
 *
 */
class SolverOperatorCreator {
public:
  virtual Solver* getSolver(const Fopr*) = 0;
  virtual Solver* getSolver(const Fopr_Herm*) = 0;//for CG like inverters
};

/////////////////////////////////////////////////
/*!
 @brief Concrete class for creating Conjugate Gradient Solver
 operator
 
 */
class SolverCGCreator : public SolverOperatorCreator {
  const XML::node Solver_node;

public:
  SolverCGCreator(const XML::node node):
    Solver_node(node){};

  Solver_CG* getSolver(const Fopr_Herm* HermitianOperator){
    return new Solver_CG(Solver_node, HermitianOperator);
  };

  Solver_CG* getSolver(const Fopr* LinearOperator){
    std::cerr<< "getSolver Error: Solver_CG requires Hermitian Operator" 
	     << std::endl;
    abort();
  };

};

/*!
 @brief Concrete class for creating BiConjugate Gradient Stabilized Solver
 operator
 
 */
class SolverBiCGStabCreator : public SolverOperatorCreator {
  const XML::node Solver_node;

public:
  SolverBiCGStabCreator(const XML::node node):
    Solver_node(node){};

  Solver_BiCGStab* getSolver(const Fopr_Herm* HermitianOperator){
    return new Solver_BiCGStab(Solver_node, HermitianOperator);
  };

  Solver_BiCGStab* getSolver(const Fopr* LinearOperator){
    return new Solver_BiCGStab(Solver_node, LinearOperator);
  };

};

///////////////////////////////////////////////////

namespace SolverOperators{
  SolverOperatorCreator* createSolverOperatorFactory(const XML::node);
}


#endif 
