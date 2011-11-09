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
#include "dirac_optimalDomainWall_4D.hpp"
#include "dirac_optimalDomainWall.hpp"
#include "Solver/solver_Factory.hpp"

/*!
 * @brief Abstract base class for creating Dirac operators
 *
 */
class DiracOperatorCreator {
public:
  virtual Dirac* getDiracOperator(Field* const) = 0;
};


/////////////////////////////////////////////////////////
/*!
 * @brief Concrete class for creating Dirac Wilson operators
 *
 */
class DiracWilsonCreator : public DiracOperatorCreator {
  const XML::node Dirac_node;

public:
  DiracWilsonCreator(const XML::node node):
    Dirac_node(node){};

  Dirac_Wilson* getDiracOperator(Field* const GaugeField){
    return new Dirac_Wilson(Dirac_node,GaugeField);
  };



};

/*!
 * @brief Concrete class for creating Dirac Optimal DWF-5d operators
 *
 */
class DiracDWF5dCreator : public DiracOperatorCreator {
  RaiiFactoryObj<DiracWilsonCreator> DiracObj;
  RaiiFactoryObj<Dirac_Wilson> Kernel;

  const XML::node Dirac_node;

public:
  DiracDWF5dCreator(XML::node node):
    Dirac_node(node){
    XML::descend(node, "Kernel");
    DiracObj.save(new DiracWilsonCreator(node)); 
    
  };

  Dirac_optimalDomainWall* getDiracOperator(Field* const GaugeField){
    Kernel.save(DiracObj.get()->getDiracOperator(GaugeField));

    return new Dirac_optimalDomainWall(Dirac_node,
				       Kernel.get());
  };

  ~DiracDWF5dCreator(){};
};


/*!
 * @brief Concrete class for creating Dirac Optimal DWF-4d operators
 *
 */
class DiracDWF4dCreator : public DiracOperatorCreator {
  RaiiFactoryObj<DiracDWF5dCreator> DiracObj;
  RaiiFactoryObj<SolverOperatorCreator> SolverDWF_Obj;

  RaiiFactoryObj<Dirac_optimalDomainWall> DWF5D_Kernel;
  RaiiFactoryObj<Fopr_DdagD> OprODWF;
  RaiiFactoryObj<Fopr_DdagD> OprPV;
  RaiiFactoryObj<Solver> SolverODWF;
  RaiiFactoryObj<Solver> SolverPV;
  
  const XML::node Dirac_node;

public:
  DiracDWF4dCreator(XML::node node):
    Dirac_node(node){
    XML::descend(node,"Kernel5d");
    DiracObj.save(new DiracDWF5dCreator(node));
    XML::next_sibling(node, "SolverDWF");
    SolverDWF_Obj.save(SolverOperators::createSolverOperatorFactory(node));
  };

  Dirac_optimalDomainWall_4D* getDiracOperator(Field* const GaugeField){
    DWF5D_Kernel.save(DiracObj.get()->getDiracOperator(GaugeField));
    OprODWF.save(new Fopr_DdagD(DWF5D_Kernel.get()));
    SolverODWF.save(SolverDWF_Obj.get()->getSolver(OprODWF.get()));
    Dirac_optimalDomainWall DPV(*DWF5D_Kernel.get(), PauliVillars);
    OprPV.save(new Fopr_DdagD(&DPV));
    SolverPV.save(SolverDWF_Obj.get()->getSolver(OprPV.get()));

    return new Dirac_optimalDomainWall_4D(*DWF5D_Kernel.get(),
					  SolverODWF.get(),
					  SolverPV.get());
  };

  ~DiracDWF4dCreator(){};
};


//Add new factories here
//....



//////////////////////////////////////////////////////////////
namespace DiracOperators {
  DiracOperatorCreator* createDiracOperatorFactory(const XML::node);

}

#endif
