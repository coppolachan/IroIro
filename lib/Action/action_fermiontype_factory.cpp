/*!
 * @file action_fermiontype_factory.cpp 
 * @brief Definition of methods for Fermion-type action factories
 * Time-stamp: <2013-08-08 10:18:36 neo>
 */
#include "action_fermiontype_factory.hpp"
#include "inputConfig.hpp"

///////////////////////////////////////////////////////////////////////////////
TwoFlavorActionFactory::TwoFlavorActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel",MANDATORY);
  DiracObj.save(Diracs::createDiracWilsonLikeFactory(node));
  XML::next_sibling(node,"Solver",MANDATORY);
  SolverObj.save(Solvers::createSolverFactory(node));
}

Action_Nf2* TwoFlavorActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  Kernel.save(DiracObj.get()->getDirac(input));
  HermitianOp.save(new Fopr_DdagD(Kernel.get()));
  Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
  
  return new Action_Nf2(input.gconf,Kernel.get(),Solv.get(),smearing, SC);
}
///////////////////////////////////////////////////////////////////////////////
NfFlavorsActionFactory::NfFlavorsActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel",MANDATORY);
  DiracObj.save(Diracs::createDiracWilsonLikeFactory(node));  
  XML::next_sibling(node,"RationalSolver",MANDATORY);
  SolverObj.save(Solvers::createRationalSolverFactory(node));
}

Action_Nf* NfFlavorsActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC) {
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));

  Kernel.save(DiracObj.get()->getDirac(input));
  HermitianOp.save(new Fopr_DdagD(Kernel.get()));
  Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
  return new Action_Nf(input.gconf,Kernel.get(),Solv.get(), 
		       Action_Nf_params(Action_node),smearing,SC);
}
///////////////////////////////////////////////////////////////////////////////
TwoFlavorRatioActionFactory::TwoFlavorRatioActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Numerator",MANDATORY);
  DiracNumObj.save(Diracs::createDiracWilsonLikeFactory(node)); 
  XML::next_sibling(node,"Denominator",MANDATORY);
  DiracDenomObj.save(Diracs::createDiracWilsonLikeFactory(node));
  XML::next_sibling(node,"SolverNumerator",MANDATORY);
  SolverNumObj.save(Solvers::createSolverFactory(node));
  XML::next_sibling(node,"SolverDenominator",MANDATORY);
  SolverDenomObj.save(Solvers::createSolverFactory(node));
}
Action_Nf2_ratio* TwoFlavorRatioActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
 InputConfig input(SC->select_conf(smearing));
  
  DiracNumerator.save(DiracNumObj.get()->getDirac(input));
  DiracDenominator.save(DiracDenomObj.get()->getDirac(input));
  
  Solver1.save(SolverNumObj.get()->getSolver(new Fopr_DdagD(DiracNumerator.get())));
  Solver2.save(SolverDenomObj.get()->getSolver(new Fopr_DdagD(DiracDenominator.get())));
  
  return new Action_Nf2_ratio(input.gconf,
			      DiracNumerator.get(),DiracDenominator.get(),
			      Solver1.get(),Solver2.get(),
			      "TwoFlavorsRatio",
			      smearing, SC); 
}
///////////////////////////////////////////////////////////////////////////////
NfFlavorRatioActionFactory::NfFlavorRatioActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Numerator",MANDATORY);
  DiracNumObj.save(Diracs::createDiracWilsonLikeFactory(node)); 
  XML::next_sibling(node,"Denominator",MANDATORY);
  DiracDenomObj.save(Diracs::createDiracWilsonLikeFactory(node));
  XML::next_sibling(node,"RationalSolverNumerator",MANDATORY);
  SolverNumObj.save(Solvers::createRationalSolverFactory(node));
  XML::next_sibling(node,"RationalSolverDenominator",MANDATORY);
  SolverDenomObj.save(Solvers::createRationalSolverFactory(node));
}
Action_Nf_ratio* NfFlavorRatioActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  DiracNumerator.save(    DiracNumObj.get()->getDirac(input));
  DiracDenominator.save(DiracDenomObj.get()->getDirac(input));
  
  Solver1.save(SolverNumObj.get()->getSolver(new Fopr_DdagD(DiracNumerator.get())));
  Solver2.save(SolverDenomObj.get()->getSolver(new Fopr_DdagD(DiracDenominator.get())));
  
  return new Action_Nf_ratio(input.gconf,
			     DiracNumerator.get(),DiracDenominator.get(),
			     Solver1.get(),Solver2.get(),
			     Action_Nf_ratio_params(Action_node),
			     "NfFlavorsRatio",
			     smearing, SC);  
}
///////////////////////////////////////////////////////////////////////////////
TwoFlavorDomainWall5dActionFactory::
TwoFlavorDomainWall5dActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel5D",MANDATORY);
  DiracObj.save(Diracs::createDiracDWF5dFactory(node));
  XML::next_sibling(node,"Solver",MANDATORY);
  SolverObj.save(Solvers::createSolverFactory(node));
}
Action_Nf2_ratio* TwoFlavorDomainWall5dActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));

  DWF5d_Kernel.save(  DiracObj.get()->getDirac(input));
  DWF5d_KernelPV.save(DiracObj.get()->getDiracPV(input));
  
  HermitianOp.save(  new Fopr_DdagD(DWF5d_Kernel.get()));
  HermitianOpPV.save(new Fopr_DdagD(DWF5d_KernelPV.get()));
  Solv.save(  SolverObj.get()->getSolver(HermitianOp.get()));
  SolvPV.save(SolverObj.get()->getSolver(HermitianOpPV.get()));
  return new Action_Nf2_ratio(input.gconf,
			      DWF5d_Kernel.get(),DWF5d_KernelPV.get(),
			      Solv.get(),SolvPV.get(),
			      "TwoFlavorsDomainWall_5D",
			      smearing, SC);
}
///////////////////////////////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
TwoFlavorDomainWall5dEO_BGQ_ActionFactory::
TwoFlavorDomainWall5dEO_BGQ_ActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node, "smeared", smearing);
  
  XML::descend(node,"Kernel5D", MANDATORY);
  DiracObj.save(new DiracEvenOdd_DWF5dFactory(node));
  XML::next_sibling(node,"Solver_DWF-EO_BGQ", MANDATORY);
  SolverObj.save(new SolverCG_DWF_opt_Factory(node));
}

Action_Nf2_ratio* TwoFlavorDomainWall5dEO_BGQ_ActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));

  DWF5d_Kernel.save(  DiracObj.get()->getDirac(input));
  DWF5d_KernelPV.save(DiracObj.get()->getDiracPV(input));
  
  Solv.save(  SolverObj.get()->getSolver(DWF5d_Kernel.get()));
  SolvPV.save(SolverObj.get()->getSolver(DWF5d_KernelPV.get()));
  
  return new Action_Nf2_ratio(input.gconf,
			      DWF5d_Kernel.get(),DWF5d_KernelPV.get(),
			      Solv.get(),SolvPV.get(),
			      "TwoFlavorsDomainWall_5D-EO_BGQ",
			      smearing, SC);
}
///////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LIBBFM
TwoFlavorDomainWall5dEO_BFM_ActionFactory::
TwoFlavorDomainWall5dEO_BFM_ActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node, "smeared", smearing);
  
  XML::descend(node,"Kernel", MANDATORY);
  DiracObj.save(new DiracBFMoperatorFactory(node));
  XML::next_sibling(node,"Solver_DWF-EO_BGQ", MANDATORY);
  SolverObj.save(new SolverCG_DWF_opt_Factory(node));
}

Action_Nf2_ratio* TwoFlavorDomainWall5dEO_BFM_ActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));

  BFM_Kernel.save( DiracObj.get()->getDirac(input));
  BFM_KernelPV.save( DiracObj.get()->getDiracPV(input));
   
  Solv.save(  SolverObj.get()->getSolver(BFM_Kernel.get()));
  SolvPV.save(SolverObj.get()->getSolver(BFM_KernelPV.get()));
  
  XML::node SolverNode = Action_node;
  XML::descend (SolverNode, "Solver_DWF-EO_BGQ",MANDATORY);  
  BFM_Kernel.get()->set_SolverParams(SolverNode);
  BFM_KernelPV.get()->set_SolverParams(SolverNode);
  
  BFM_Kernel.get()->initialize();
  BFM_KernelPV.get()->initialize();

  return new Action_Nf2_ratio(input.gconf,
			      BFM_Kernel.get(),BFM_KernelPV.get(),
			      Solv.get(),SolvPV.get(),
			      "TwoFlavorsDomainWall_5D-EO_BFM",
			      smearing, SC);
  
}
#endif
///////////////////////////////////////////////////////////////////////////////
TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory::
TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);

  XML::descend(node,"Numerator", MANDATORY);
  DiracNumObj.save(new DiracEvenOdd_DWF5dFactory(node));
  XML::next_sibling(node,"Denominator", MANDATORY);
  DiracDenomObj.save(new DiracEvenOdd_DWF5dFactory(node));
  XML::next_sibling(node,"SolverNumerator", MANDATORY);
  SolverNumObj.save(new SolverCG_DWF_opt_Factory(node));
  XML::next_sibling(node,"SolverDenominator", MANDATORY);
  SolverDenomObj.save(new SolverCG_DWF_opt_Factory(node));
}
Action_Nf2_ratio* TwoFlavorRatioDomainWall5dEO_BGQ_ActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  DiracNumerator.save(DiracNumObj.get()->getDirac(input));
  DiracDenominator.save(DiracDenomObj.get()->getDirac(input));
  
  Solver1.save(SolverNumObj.get()->getSolver(DiracNumerator.get()));
  Solver2.save(SolverDenomObj.get()->getSolver(DiracDenominator.get()));
  
  return new Action_Nf2_ratio(input.gconf,
			      DiracNumerator.get(),DiracDenominator.get(),
			      Solver1.get(),Solver2.get(),
			      "TwoFlavorsRatioDomainWall_5D-EO_BGQ",
			      smearing, SC); 
}
#endif
///////////////////////////////////////////////////////////////////////////////
NfFlavorDomainWall5dActionFactory::
NfFlavorDomainWall5dActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel5D",MANDATORY);
  DiracObj.save(Diracs::createDiracDWF5dFactory(node));
  XML::next_sibling(node,"RationalSolver",MANDATORY);
  SolverObj.save(Solvers::createRationalSolverFactory(node));
}
Action_Nf_ratio* NfFlavorDomainWall5dActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  DWF5d_Kernel.save(  DiracObj.get()->getDirac(input));
  DWF5d_KernelPV.save(DiracObj.get()->getDiracPV(input));
  
  HermitianOp.save(  new Fopr_DdagD(DWF5d_Kernel.get()));
  HermitianOpPV.save(new Fopr_DdagD(DWF5d_KernelPV.get()));
  Solv.save(  SolverObj.get()->getSolver(HermitianOp.get()));
  SolvPV.save(SolverObj.get()->getSolver(HermitianOpPV.get()));
  return new Action_Nf_ratio(input.gconf,
			     DWF5d_Kernel.get(),DWF5d_KernelPV.get(),
			     Solv.get(),SolvPV.get(),
			     Action_Nf_ratio_params(Action_node),
			     "NfFlavorsDomainWall_5D",
			     smearing,SC);
}
///////////////////////////////////////////////////////////////////////////////
#ifdef IBM_BGQ_WILSON
NfFlavorDomainWall5d_EO_BGQ_ActionFactory::
NfFlavorDomainWall5d_EO_BGQ_ActionFactory(XML::node node)
:Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel5D",MANDATORY);
  DiracObj.save(new DiracEvenOdd_DWF5dFactory(node));
  XML::next_sibling(node,"RationalSolver",MANDATORY);
  SolverObj.save(new RationalSolverCGFactory_DWF_Optimized(node));
}

Action_Nf_ratio* NfFlavorDomainWall5d_EO_BGQ_ActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  DWF5d_Kernel.save(  DiracObj.get()->getDirac(input));
  DWF5d_KernelPV.save(DiracObj.get()->getDiracPV(input));
  
  Solv.save(  SolverObj.get()->getSolver(DWF5d_Kernel.get()));
  SolvPV.save(SolverObj.get()->getSolver(DWF5d_KernelPV.get()));
  return new Action_Nf_ratio(input.gconf,
			     DWF5d_Kernel.get(),DWF5d_KernelPV.get(),
			     Solv.get(),SolvPV.get(),
			     Action_Nf_ratio_params(Action_node),
			     "NfFlavorsDomainWall_5D-EO_BGQ",
			     smearing,SC);
}
#endif
///////////////////////////////////////////////////////////////////////////////
FourFlavorStaggeredActionFactory::
FourFlavorStaggeredActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Kernel",MANDATORY);
  DiracObj.save(Diracs::createDiracStaggeredEvenOddLikeFactory(node));
  XML::next_sibling(node,"Solver",MANDATORY);
  SolverObj.save(Solvers::createSolverFactory(node));
}

Action_staggered* FourFlavorStaggeredActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  Kernel.save(DiracObj.get()->getDirac(input));
  HermitianOp.save(new Fopr_HD(Kernel.get()));
  Solv.save(SolverObj.get()->getSolver(HermitianOp.get()));
  
  return new Action_staggered(input.gconf,Kernel.get(),Solv.get(),smearing,SC);
}
///////////////////////////////////////////////////////////////////////////////
FourFlavorStaggeredRatioActionFactory::
FourFlavorStaggeredRatioActionFactory(XML::node node)
  :Action_node(node),smearing(false){
  XML::read(node,"smeared",smearing);
  
  XML::descend(node,"Numerator",MANDATORY);
  DiracNumObj.save(Diracs::createDiracStaggeredEvenOddLikeFactory(node)); 
  XML::next_sibling(node,"Denominator",MANDATORY);
  DiracDenomObj.save(Diracs::createDiracStaggeredEvenOddLikeFactory(node));
  XML::next_sibling(node,"SolverNumerator",MANDATORY);
  SolverNumObj.save(Solvers::createSolverFactory(node));
  XML::next_sibling(node,"SolverDenominator",MANDATORY);
  SolverDenomObj.save(Solvers::createSolverFactory(node));
}

Action_staggered_ratio* FourFlavorStaggeredRatioActionFactory::
getFermionAction(GaugeField* const F,SmartConf* const SC){
  // select links according to smearing
  InputConfig input(SC->select_conf(smearing));
  
  DiracNumerator_ee.save(DiracNumObj.get()->getDirac(input));
  DiracNumerator_oo.save(DiracNumObj.get()->getDoo(input));
  DiracDenominator.save(DiracDenomObj.get()->getDirac(input));
  
  Solver1e.save(SolverNumObj.get()->getSolver(new Fopr_HD(DiracNumerator_ee.get())));
  Solver1o.save(SolverNumObj.get()->getSolver(new Fopr_HD(DiracNumerator_oo.get())));
  Solver2e.save(SolverDenomObj.get()->getSolver(new Fopr_HD(DiracDenominator.get())));
  
  return new Action_staggered_ratio(input.gconf,
				    DiracNumerator_ee.get(),DiracDenominator.get(),
				    Solver1e.get(),Solver1o.get(),Solver2e.get(),
				    smearing,SC); 
}
