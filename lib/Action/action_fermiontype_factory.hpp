/*!
 * @file action_fermiontype_factory.hpp 
 *
 * @brief Declaration of FermionType action factories
 */
#ifndef ACTION_FERMION_FACT_
#define ACTION_FERMION_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

#include "action_fermiontype_factory_abs.hpp"
#include "Action/action_Nf2.hpp"
#include "Action/action_Nf.hpp"
#include "Action/action_Nf2_ratio.hpp"
#include "Action/action_Nf_ratio.hpp"
#include "Action/action_Nf2_DomainWall.hpp"
#include "Solver/solver_CG.hpp"
#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Solver/solver_Factory.hpp"


///////////////////////////////////////////////////////////////////////

class TwoFlavorActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Fopr_DdagD> HermitianOp;
  RaiiFactoryObj<Solver> Solv;
 
  const XML::node Action_node;

public:
  TwoFlavorActionFactory(XML::node node):Action_node(node){
    XML::descend(node,"Kernel",MANDATORY);
    DiracObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));
    XML::next_sibling(node,"Solver", MANDATORY);
    SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
  }

  ~TwoFlavorActionFactory(){}
private:
  Action_Nf2* getFermionAction(GaugeField* const F){
    Kernel.save(DiracObj.get()->getDiracOperator(&(F->data)));
    HermitianOp.save(new Fopr_DdagD(Kernel.get()));
    Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
    return new Action_Nf2(F, Kernel.get(), Solv.get());
  }
};

//////////////////////////////////////////////////////////////////

class NfFlavorsActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracObj;
  RaiiFactoryObj<RationalSolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Fopr_DdagD> HermitianOp;
  RaiiFactoryObj<RationalSolver> Solv;

  const XML::node Action_node;
  
public:
  NfFlavorsActionFactory(XML::node node):Action_node(node) {
    XML::descend(node, "Kernel", MANDATORY);
    DiracObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));  
    XML::next_sibling(node,"RationalSolver", MANDATORY);
    SolverObj.save(SolverOperators::createRationalSolverOperatorFactory(node));
  }

  ~NfFlavorsActionFactory(){}
private: 
  Action_Nf* getFermionAction(GaugeField* const F) {
    Kernel.save(DiracObj.get()->getDiracOperator(&(F->data)));
    HermitianOp.save(new Fopr_DdagD(Kernel.get()));
    Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
    return new Action_Nf(F, Kernel.get(), Solv.get(), 
			 Action_Nf_params(Action_node));
  }

};

////////////////////////////////////////////////////

class TwoFlavorRatioActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracNumObj;
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracDenomObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverNumObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverDenomObj;

  RaiiFactoryObj<DiracWilsonLike> DiracNumerator;
  RaiiFactoryObj<DiracWilsonLike> DiracDenominator;
  RaiiFactoryObj<Solver> Solver1;
  RaiiFactoryObj<Solver> Solver2;

  const XML::node Action_node;

public:
  ~TwoFlavorRatioActionFactory(){}

  TwoFlavorRatioActionFactory(XML::node node):Action_node(node){
    XML::descend(node,"Numerator", MANDATORY);
    DiracNumObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node)); 
    XML::next_sibling(node,"Denominator", MANDATORY);
    DiracDenomObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));
    XML::next_sibling(node,"SolverNumerator", MANDATORY);
    SolverNumObj.save(SolverOperators::createSolverOperatorFactory(node));
    XML::next_sibling(node,"SolverDenominator", MANDATORY);
    SolverDenomObj.save(SolverOperators::createSolverOperatorFactory(node));
  }
  
private:  
  Action_Nf2_ratio* getFermionAction(GaugeField* const F){
    DiracNumerator.save(DiracNumObj.get()->getDiracOperator(&(F->data)));
    DiracDenominator.save(DiracDenomObj.get()->getDiracOperator(&(F->data)));
    
    Solver1.save(SolverNumObj.get()->getSolver(new Fopr_DdagD(DiracNumerator.get())));
    Solver2.save(SolverDenomObj.get()->getSolver(new Fopr_DdagD(DiracDenominator.get())));

    return new Action_Nf2_ratio(F,
				DiracNumerator.get(),
				DiracDenominator.get(),
				Solver1.get(),
				Solver2.get()); 
  }
};

////////////////////////////////////////////////////

class NfFlavorRatioActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracNumObj;
  RaiiFactoryObj<DiracWilsonLikeOperatorFactory> DiracDenomObj;
  RaiiFactoryObj<RationalSolverOperatorFactory> SolverNumObj;
  RaiiFactoryObj<RationalSolverOperatorFactory> SolverDenomObj;

  RaiiFactoryObj<DiracWilsonLike> DiracNumerator;
  RaiiFactoryObj<DiracWilsonLike> DiracDenominator;
  RaiiFactoryObj<RationalSolver> Solver1;
  RaiiFactoryObj<RationalSolver> Solver2;

  const XML::node Action_node;

public:
  ~NfFlavorRatioActionFactory(){}

  NfFlavorRatioActionFactory(XML::node node):Action_node(node){
    XML::descend(node,"Numerator", MANDATORY);
    DiracNumObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node)); 
    XML::next_sibling(node,"Denominator", MANDATORY);
    DiracDenomObj.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));
    XML::next_sibling(node,"RationalSolverNumerator", MANDATORY);
    SolverNumObj.save(SolverOperators::createRationalSolverOperatorFactory(node));
    XML::next_sibling(node,"RationalSolverDenominator", MANDATORY);
    SolverDenomObj.save(SolverOperators::createRationalSolverOperatorFactory(node));
  }
  
private:  
  Action_Nf_ratio* getFermionAction(GaugeField* const F){
    DiracNumerator.save(DiracNumObj.get()->getDiracOperator(&(F->data)));
    DiracDenominator.save(DiracDenomObj.get()->getDiracOperator(&(F->data)));
    
    Solver1.save(SolverNumObj.get()->getSolver(new Fopr_DdagD(DiracNumerator.get())));
    Solver2.save(SolverDenomObj.get()->getSolver(new Fopr_DdagD(DiracDenominator.get())));

    return new Action_Nf_ratio(F,
			       DiracNumerator.get(),
			       DiracDenominator.get(),
			       Solver1.get(),
			       Solver2.get(),
			       Action_Nf_ratio_params(Action_node)); 
  }
};



////////////////////////////////////////////////////

class TwoFlavorDomainWall5dActionFactory : public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dOperatorFactory> DiracObj;
  RaiiFactoryObj<SolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> DWF5d_Kernel;
  RaiiFactoryObj<DiracWilsonLike> DWF5d_KernelPV;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOp;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOpPV;
  RaiiFactoryObj<Solver> Solv;
  RaiiFactoryObj<Solver> SolvPV;
  
  const XML::node Action_node;

public:
  ~TwoFlavorDomainWall5dActionFactory(){}

  TwoFlavorDomainWall5dActionFactory(XML::node node):Action_node(node){
    XML::descend(node,"Kernel5D", MANDATORY);
    DiracObj.save(DiracOperators::createDiracDWF5dOperatorFactory(node));
    XML::next_sibling(node,"Solver", MANDATORY);
    SolverObj.save(SolverOperators::createSolverOperatorFactory(node));
  }

private:  
  Action_Nf2_ratio* getFermionAction(GaugeField* const F){
    DWF5d_Kernel.save(  DiracObj.get()->getDiracOperator(&(F->data)));
    DWF5d_KernelPV.save(DiracObj.get()->getDiracOperatorPV(&(F->data)));

    HermitianOp.save(  new Fopr_DdagD_Precondition(DWF5d_Kernel.get()));
    HermitianOpPV.save(new Fopr_DdagD_Precondition(DWF5d_KernelPV.get()));
    Solv.save(  SolverObj.get()->getSolver(HermitianOp.get()));
    SolvPV.save(SolverObj.get()->getSolver(HermitianOpPV.get()));
    return new Action_Nf2_ratio(F,
				DWF5d_Kernel.get(),
				DWF5d_KernelPV.get(),
				Solv.get(),
				SolvPV.get());
  }
};

////////////////////////////////////////////////////

class NfFlavorDomainWall5dActionFactory : public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dOperatorFactory> DiracObj;
  RaiiFactoryObj<RationalSolverOperatorFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> DWF5d_Kernel;
  RaiiFactoryObj<DiracWilsonLike> DWF5d_KernelPV;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOp;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOpPV;
  RaiiFactoryObj<RationalSolver> Solv;
  RaiiFactoryObj<RationalSolver> SolvPV;
  
  const XML::node Action_node;

public:
  ~NfFlavorDomainWall5dActionFactory(){}

  NfFlavorDomainWall5dActionFactory(XML::node node):Action_node(node){
    XML::descend(node,"Kernel5D", MANDATORY);
    DiracObj.save(DiracOperators::createDiracDWF5dOperatorFactory(node));
    XML::next_sibling(node,"RationalSolver", MANDATORY);
    SolverObj.save(SolverOperators::createRationalSolverOperatorFactory(node));
  }

private:  
  Action_Nf_ratio* getFermionAction(GaugeField* const F){
    DWF5d_Kernel.save(  DiracObj.get()->getDiracOperator(&(F->data)));
    DWF5d_KernelPV.save(DiracObj.get()->getDiracOperatorPV(&(F->data)));

    HermitianOp.save(  new Fopr_DdagD_Precondition(DWF5d_Kernel.get()));
    HermitianOpPV.save(new Fopr_DdagD_Precondition(DWF5d_KernelPV.get()));
    Solv.save(  SolverObj.get()->getSolver(HermitianOp.get()));
    SolvPV.save(SolverObj.get()->getSolver(HermitianOpPV.get()));
    return new Action_Nf_ratio(F,
			       DWF5d_Kernel.get(),
			       DWF5d_KernelPV.get(),
			       Solv.get(),
			       SolvPV.get(),
			       Action_Nf_ratio_params(Action_node));
  }
};


//Add new factories here
//....

#endif
