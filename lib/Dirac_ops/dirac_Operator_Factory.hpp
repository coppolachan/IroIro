/*!
 * @file dirac_Factory.hpp 
 *
 * @brief Declaration of Dirac operators factories
 */

#ifndef DIRAC_FACT_
#define DIRAC_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"
#include "dirac_wilson.h"
#include "dirac_wilson_EvenOdd.h"
#include "dirac_clover.hpp"
#include "dirac_DomainWall_4D.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "Solver/solver_Factory.hpp"

/*!
 * @brief Abstract base class for creating Dirac operators
 *
 */
class DiracWilsonLikeOperatorFactory {
public:
  virtual DiracWilsonLike* getDiracOperator(Field* const) = 0;
};

/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Wilson operators
 *
 */
class DiracWilsonFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;

public:
  DiracWilsonFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson* getDiracOperator(Field* const GaugeField){
    return new Dirac_Wilson(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Wilson EvenOdd operators
 *
 */
class DiracWilsonEvenOddFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;

public:
  DiracWilsonEvenOddFactory(const XML::node node):Dirac_node(node){}

  Dirac_Wilson_EvenOdd* getDiracOperator(Field* const GaugeField){
    return new Dirac_Wilson_EvenOdd(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Clover operators
 *
 */
class DiracCloverFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;

public:
  DiracCloverFactory(const XML::node node):Dirac_node(node){}

  Dirac_Clover* getDiracOperator(Field* const GaugeField){
    return new Dirac_Clover(Dirac_node,GaugeField);
  }
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Optimal DWF-5d operators
 */
class DiracDWF5dFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall* getDiracOperator(Field* const GaugeField){
    return new Dirac_optimalDomainWall(Dirac_node,GaugeField);
  }
  ~DiracDWF5dFactory(){}
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Optimal DWF-5d e/o operators
 */
class DiracDWF5dEvenOddFactory : public DiracWilsonLikeOperatorFactory {
  const XML::node Dirac_node;
public:
  DiracDWF5dEvenOddFactory(XML::node node):Dirac_node(node){}

  Dirac_optimalDomainWall_EvenOdd* getDiracOperator(Field* const GaugeField){
    return new Dirac_optimalDomainWall_EvenOdd(Dirac_node,GaugeField);
  }
  ~DiracDWF5dEvenOddFactory(){}
};
/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Optimal DWF-4d operators
 */
class DiracDWF4dFactory : public DiracWilsonLikeOperatorFactory {
  RaiiFactoryObj<DiracDWF5dFactory> DiracObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverDWF_Obj;

  RaiiFactoryObj<Dirac_optimalDomainWall> DWF5D_Kernel;
  RaiiFactoryObj<Dirac_optimalDomainWall> DPV;
  RaiiFactoryObj<Fopr_DdagD_Precondition> OprODWF;
  RaiiFactoryObj<Fopr_DdagD_Precondition> OprPV;
  RaiiFactoryObj<Solver> SolverODWF;
  RaiiFactoryObj<Solver> SolverPV;
  
  const XML::node Dirac_node;

public:
  DiracDWF4dFactory(XML::node node)
    :Dirac_node(node){
    XML::descend(node,"Kernel5d", MANDATORY);
    DiracObj.save(new DiracDWF5dFactory(node));
    XML::next_sibling(node, "SolverDWF", MANDATORY);
    SolverDWF_Obj.save(SolverOperators::createSolverOperatorFactory(node));
  }

  Dirac_optimalDomainWall_4D* getDiracOperator(Field* const GaugeField){
    DWF5D_Kernel.save(DiracObj.get()->getDiracOperator(GaugeField));
    OprODWF.save(new Fopr_DdagD_Precondition(DWF5D_Kernel.get()));
    SolverODWF.save(SolverDWF_Obj.get()->getSolver(OprODWF.get()));
    DPV.save(new Dirac_optimalDomainWall(*DWF5D_Kernel.get(), 
					 DomainWallFermions::PauliVillars_tag()));
    OprPV.save(new Fopr_DdagD_Precondition(DPV.get()));
    SolverPV.save(SolverDWF_Obj.get()->getSolver(OprPV.get()));

    return new Dirac_optimalDomainWall_4D(*DWF5D_Kernel.get(),
					  SolverODWF.get(),
					  SolverPV.get());
  }
  ~DiracDWF4dFactory(){}
};


//Add new factories here
//....

//////////////////////////////////////////////////////////////
namespace DiracOperators {
  DiracWilsonLikeOperatorFactory* 
  createDiracWilsonLikeOperatorFactory(const XML::node);

}

#endif
