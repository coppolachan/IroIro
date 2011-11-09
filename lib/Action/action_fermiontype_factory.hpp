/*!
 * @file action_fermiontype_factory.hpp 
 *
 * @brief Declaration of FermionType action factories
 */

#ifndef ACTION_FERMION_FACT_
#define ACTION_FERMION_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

#include "action_Factory.h"
#include "Action/action_Nf2.h"
#include "Action/action_Nf2_ratio.h"
#include "Solver/solver_CG.h"
#include "Dirac_ops/dirac_Factory.hpp"
#include "Solver/solver_Factory.hpp"

class FermionActionCreator : public ActionCreator {
  virtual Action* getFermionAction(const Format::Format_G&,
				   Field* const) = 0;
public:
  Action* getAction(const Format::Format_G& GaugeForm,
		    Field* const GaugeField) {
    return getFermionAction(GaugeForm,
			    GaugeField);
  }

};

///////////////////////////////////////////////////////////////////////

class TwoFlavorActionCreator : public FermionActionCreator {
  RaiiFactoryObj<DiracOperatorCreator> DiracObj;
  RaiiFactoryObj<SolverOperatorCreator> SolverObj;

  RaiiFactoryObj<Dirac> Kernel;
  RaiiFactoryObj<Fopr_DdagD> HermitianOp;
  RaiiFactoryObj<Solver> Solv;
 
  const XML::node Action_node;

public:
  TwoFlavorActionCreator(XML::node node):
    Action_node(node)
    {
      XML::descend(node,"Kernel");
      DiracObj.save(DiracOperators::createDiracOperatorFactory(node));
      XML::next_sibling(node,"Solver");
      SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
    };

  ~TwoFlavorActionCreator(){}
private:

  Action_Nf2* getFermionAction(const Format::Format_G& Form,
			       Field* const GaugeField){
    Kernel.save(DiracObj.get()->getDiracOperator(GaugeField));
    HermitianOp.save(new Fopr_DdagD(Kernel.get()));
    Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
    return new Action_Nf2(GaugeField, Kernel.get(), Solv.get());

  };
  
};

class TwoFlavorRatioActionCreator : public FermionActionCreator {
  RaiiFactoryObj<DiracOperatorCreator> DiracNumObj;
  RaiiFactoryObj<DiracOperatorCreator> DiracDenomObj;
  RaiiFactoryObj<SolverOperatorCreator> SolverNumObj;
  RaiiFactoryObj<SolverOperatorCreator> SolverDenomObj;

  RaiiFactoryObj<Dirac> DiracNumerator;
  RaiiFactoryObj<Dirac> DiracDenominator;
  RaiiFactoryObj<Solver> Solver1;
  RaiiFactoryObj<Solver> Solver2;

  const XML::node Action_node;

public:
  ~TwoFlavorRatioActionCreator(){}

  TwoFlavorRatioActionCreator(XML::node node):
    Action_node(node){
    XML::descend(node,"Numerator");
    DiracNumObj.save(DiracOperators::createDiracOperatorFactory(node)); 
    XML::next_sibling(node,"Denominator");
    DiracDenomObj.save(DiracOperators::createDiracOperatorFactory(node));
    XML::next_sibling(node,"SolverNumerator");
    SolverNumObj.save(SolverOperators::createSolverOperatorFactory(node));
    XML::next_sibling(node,"SolverDenominator");
    SolverDenomObj.save(SolverOperators::createSolverOperatorFactory(node));

};


  
private:  
  Action_Nf2_ratio* getFermionAction(const Format::Format_G& Form,
				     Field* const GaugeField){
    DiracNumerator.save(DiracNumObj.get()->getDiracOperator(GaugeField));
    DiracDenominator.save(DiracDenomObj.get()->getDiracOperator(GaugeField));
    
    Solver1.save(SolverNumObj.get()->getSolver(new Fopr_DdagD(DiracNumerator.get())));
    Solver2.save(SolverDenomObj.get()->getSolver(new Fopr_DdagD(DiracDenominator.get())));

    return new Action_Nf2_ratio(GaugeField,
				DiracNumerator.get(),
				DiracDenominator.get(),
				Solver1.get(),
				Solver2.get()); 
  };
  

};

//Add new factories here
//....


////////////////////////////////////////////////////
namespace FermionAction {
  FermionActionCreator* createFermionActionFactory(XML::node);

}

#endif
