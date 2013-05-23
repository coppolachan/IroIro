/*!
 * @file action_fermiontype_factory.hpp 
 * @brief Declaration of Fermion-type action factories
 * Time-stamp: <2013-05-23 11:16:08 noaki>
 */
#ifndef ACTION_FERMION_FACT_
#define ACTION_FERMION_FACT_

#include "include/pugi_interface.h"
#include "Tools/RAIIFactory.hpp"

// Include here the files containing 
// the declarations of the actions
#include "action_fermiontype_factory_creator.hpp"
#include "Action/action_Nf2.hpp"
#include "Action/action_Nf.hpp"
#include "Action/action_Nf2_ratio.hpp"
#include "Action/action_Nf_ratio.hpp"
#include "Action/action_staggered.hpp"
#include "Action/action_staggered_ratio.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"

///////////////////////////////////////////////////////////////////////
class TwoFlavorActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracObj;
  RaiiFactoryObj<SolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Fopr_DdagD> HermitianOp;
  RaiiFactoryObj<Solver> Solv;
 
  const XML::node Action_node;
  bool smearing;

  Action_Nf2* getFermionAction(GaugeField* const, 
			       SmartConf* const);
public:
  TwoFlavorActionFactory(XML::node);
  ~TwoFlavorActionFactory(){}
};

//////////////////////////////////////////////////////////////////
class NfFlavorsActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracObj;
  RaiiFactoryObj<RationalSolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> Kernel;
  RaiiFactoryObj<Fopr_DdagD> HermitianOp;
  RaiiFactoryObj<RationalSolver> Solv;

  const XML::node Action_node;
  bool smearing;
  
  Action_Nf* getFermionAction(GaugeField* const,SmartConf* const);
public:
  NfFlavorsActionFactory(XML::node);
  ~NfFlavorsActionFactory(){}
};

////////////////////////////////////////////////////
class TwoFlavorRatioActionFactory :public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracNumObj;
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracDenomObj;
  RaiiFactoryObj<SolverFactory> SolverNumObj;
  RaiiFactoryObj<SolverFactory> SolverDenomObj;

  RaiiFactoryObj<DiracWilsonLike> DiracNumerator;
  RaiiFactoryObj<DiracWilsonLike> DiracDenominator;
  RaiiFactoryObj<Solver> Solver1;
  RaiiFactoryObj<Solver> Solver2;

  const XML::node Action_node;
  bool smearing;

  Action_Nf2_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~TwoFlavorRatioActionFactory(){}
  TwoFlavorRatioActionFactory(XML::node);
};

////////////////////////////////////////////////////
class NfFlavorRatioActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracNumObj;
  RaiiFactoryObj<DiracWilsonLikeFactory> DiracDenomObj;
  RaiiFactoryObj<RationalSolverFactory> SolverNumObj;
  RaiiFactoryObj<RationalSolverFactory> SolverDenomObj;

  RaiiFactoryObj<DiracWilsonLike> DiracNumerator;
  RaiiFactoryObj<DiracWilsonLike> DiracDenominator;
  RaiiFactoryObj<RationalSolver> Solver1;
  RaiiFactoryObj<RationalSolver> Solver2;

  const XML::node Action_node;
  bool smearing;

  Action_Nf_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~NfFlavorRatioActionFactory(){}
  NfFlavorRatioActionFactory(XML::node);
};

////////////////////////////////////////////////////
class TwoFlavorDomainWall5dActionFactory :public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dFactory> DiracObj;
  RaiiFactoryObj<SolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> DWF5d_Kernel;
  RaiiFactoryObj<DiracWilsonLike> DWF5d_KernelPV;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOp;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOpPV;
  RaiiFactoryObj<Solver> Solv;
  RaiiFactoryObj<Solver> SolvPV;
  
  const XML::node Action_node;
  bool smearing;

  Action_Nf2_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~TwoFlavorDomainWall5dActionFactory(){}
  TwoFlavorDomainWall5dActionFactory(XML::node);
};

////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
class TwoFlavorDomainWall5dEO_BGQ_ActionFactory : public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dEvenOddFactory> DiracObj;
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverObj;

  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DWF5d_Kernel;
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DWF5d_KernelPV;
  RaiiFactoryObj<Solver> Solv;
  RaiiFactoryObj<Solver> SolvPV;
  
  const XML::node Action_node;
  bool smearing;

  Action_Nf2_ratio* getFermionAction(GaugeField* const,
				     SmartConf* const);
public:
  ~TwoFlavorDomainWall5dEO_BGQ_ActionFactory(){}
  TwoFlavorDomainWall5dEO_BGQ_ActionFactory(XML::node);
};

////////////////////////////////////////////////////
class TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory: public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dEvenOddFactory> DiracNumObj;
  RaiiFactoryObj<DiracDWF5dEvenOddFactory> DiracDenomObj;
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverNumObj;
  RaiiFactoryObj<SolverCG_DWF_opt_Factory> SolverDenomObj;

  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DiracNumerator;
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DiracDenominator;
  RaiiFactoryObj<Solver> Solver1;
  RaiiFactoryObj<Solver> Solver2;

  const XML::node Action_node;
  bool smearing;

  Action_Nf2_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory(){}
  TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory(XML::node);
};

#endif

////////////////////////////////////////////////////
class NfFlavorDomainWall5dActionFactory : public FermionActionFactory {

  RaiiFactoryObj<DiracDWF5dFactory> DiracObj;
  RaiiFactoryObj<RationalSolverFactory> SolverObj;

  RaiiFactoryObj<DiracWilsonLike> DWF5d_Kernel;
  RaiiFactoryObj<DiracWilsonLike> DWF5d_KernelPV;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOp;
  RaiiFactoryObj<Fopr_DdagD_Precondition> HermitianOpPV;
  RaiiFactoryObj<RationalSolver> Solv;
  RaiiFactoryObj<RationalSolver> SolvPV;
  
  const XML::node Action_node;
  bool smearing;

  Action_Nf_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~NfFlavorDomainWall5dActionFactory(){}
  NfFlavorDomainWall5dActionFactory(XML::node);
private:  

};
////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
class NfFlavorDomainWall5d_EO_BGQ_ActionFactory : public FermionActionFactory {
  RaiiFactoryObj<DiracDWF5dEvenOddFactory> DiracObj;
  RaiiFactoryObj<RationalSolverCGFactory_DWF_Optimized> SolverObj;

  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DWF5d_Kernel;
  RaiiFactoryObj<Dirac_optimalDomainWall_EvenOdd> DWF5d_KernelPV;
  RaiiFactoryObj<RationalSolver_DWF_Optimized> Solv;
  RaiiFactoryObj<RationalSolver_DWF_Optimized> SolvPV;
  
  const XML::node Action_node;
  bool smearing;

  Action_Nf_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  ~NfFlavorDomainWall5d_EO_BGQ_ActionFactory(){}
  NfFlavorDomainWall5d_EO_BGQ_ActionFactory(XML::node);
};
#endif
////////////////////////////////////////////////////
class FourFlavorStaggeredActionFactory :public FermionActionFactory{
  RaiiFactoryObj<DiracStaggeredEvenOddLikeFactory> DiracObj;
  RaiiFactoryObj<SolverFactory> SolverObj;

  RaiiFactoryObj<DiracStaggeredEvenOddLike> Kernel;
  RaiiFactoryObj<Fopr_HD> HermitianOp;
  RaiiFactoryObj<Solver> Solv;
 
  const XML::node Action_node;
  bool smearing;

  Action_staggered* getFermionAction(GaugeField* const,SmartConf* const);
public:
  FourFlavorStaggeredActionFactory(XML::node);  
  ~FourFlavorStaggeredActionFactory(){}
};
////////////////////////////////////////////////////
class FourFlavorStaggeredRatioActionFactory :public FermionActionFactory{
  RaiiFactoryObj<DiracStaggeredEvenOddLikeFactory> DiracNumObj;
  RaiiFactoryObj<DiracStaggeredEvenOddLikeFactory> DiracDenomObj;
  RaiiFactoryObj<SolverFactory> SolverNumObj;
  RaiiFactoryObj<SolverFactory> SolverDenomObj;

  RaiiFactoryObj<DiracStaggeredEvenOddLike> DiracNumerator_ee;
  RaiiFactoryObj<DiracStaggeredEvenOddLike> DiracNumerator_oo;
  RaiiFactoryObj<DiracStaggeredEvenOddLike> DiracDenominator;
  RaiiFactoryObj<Solver> Solver1e;
  RaiiFactoryObj<Solver> Solver1o;
  RaiiFactoryObj<Solver> Solver2e;

  const XML::node Action_node;
  bool smearing;

  Action_staggered_ratio* getFermionAction(GaugeField* const,SmartConf* const);
public:
  FourFlavorStaggeredRatioActionFactory(XML::node);
  ~FourFlavorStaggeredRatioActionFactory(){}
};

// Add new factories here and in the .cpp file
// ....

#endif
