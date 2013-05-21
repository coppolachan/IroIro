/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-05-21 11:36:58 noaki>
 */
#include "dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"
#include "Communicator/comm_io.hpp"
/*
#include "dirac_wilson.hpp"
#include "dirac_wilson_EvenOdd.hpp"
#include "dirac_wilson_Brillouin.hpp"
#include "dirac_clover.hpp"
#include "dirac_DomainWall_4D_fullSolv.hpp"
#include "dirac_DomainWall_4D_eoSolv.hpp"
#include "dirac_DomainWall.hpp"
#include "dirac_DomainWall_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd.hpp"
#include "dirac_staggered_EvenOdd_Adjoint.hpp"
*/
#include "eoUtils.hpp"
#include <string.h>

/// Dirac_Wilson
Dirac_Wilson* DiracWilsonFactory::getDiracOperatorWL(Field* const Gfield){
  return new Dirac_Wilson(Dirac_node_,Gfield);
}
/// Dirac_Wilson_EvenOdd
Dirac_Wilson_EvenOdd* DiracWilsonEvenOddFactory::
getDiracOperatorEO(Field* const Gfield){
  return new Dirac_Wilson_EvenOdd(Dirac_node_,Gfield);
}

/// Dirac_Wilson_Brillouin
Dirac_Wilson_Brillouin* DiracWilsonBrillouinFactory::
getDiracOperatorWL(Field* const Gfield){
  return new Dirac_Wilson_Brillouin(Dirac_node_,Gfield,type_);}

/// Dirac_Clover
Dirac_Clover* DiracCloverFactory::getDiracOperatorWL(Field* const Gfield){
  return new Dirac_Clover(Dirac_node_,Gfield); }

/// Dirac_Mobius
DiracMobiusFactory::DiracMobiusFactory(XML::node node)
  :Dirac_node_(node){
  XML::node dnode = node;
  XML::descend(dnode,"Dirac_denominator",MANDATORY);
  XML::descend(node,"Solver",MANDATORY);

  DiracFactory_.save(DiracOperators::createDiracWilsonLikeOperatorFactory(dnode));
  SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
}

Dirac_Mobius* DiracMobiusFactory::getDiracOperatorWL(Field* const Gfield){
  D_.save(DiracFactory_.get()->getDiracOperatorWL(Gfield));
  Fopr_.save(new Fopr_DdagD_Precondition(D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));
  return new Dirac_Mobius(Dirac_node_,D_.get(),Solver_.get());
}

/// Dirac_optimalDomainWall
DiracDWF5dFactory::DiracDWF5dFactory(XML::node node):Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(DiracOperators::createDiracWilsonLikeOperatorFactory(node));
}

Dirac_optimalDomainWall* DiracDWF5dFactory::getDiracOperatorWL(Field* const Gfield){
  return new Dirac_optimalDomainWall(Dirac_node_,
				     KernelFactory_.get()->getDiracOperatorWL(Gfield),
				     Gfield); 
}

Dirac_optimalDomainWall* DiracDWF5dFactory::getDiracOperatorPV(Field* const Gfield){
  return new Dirac_optimalDomainWall(Dirac_node_,
				     KernelFactory_.get()->getDiracOperatorWL(Gfield),
				     Gfield,PauliVillars); 
}

/// Dirac_optimalDomainWall_EvenOdd
DiracDWF5dEvenOddFactory::DiracDWF5dEvenOddFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(DiracOperators::createDiracWilsonLikeEvenOddOperatorFactory(node));
}

Dirac_optimalDomainWall_EvenOdd* DiracDWF5dEvenOddFactory::
getDiracOperatorWL(Field* const Gfield){
  Kernel_.save(KernelFactory_.get()->getDiracOperatorEO(Gfield));
  return new Dirac_optimalDomainWall_EvenOdd(Dirac_node_,Kernel_.get(),Gfield); 
}

Dirac_optimalDomainWall_EvenOdd* DiracDWF5dEvenOddFactory::
getDiracOperatorPV(Field* const Gfield){
  Kernel_.save(KernelFactory_.get()->getDiracOperatorEO(Gfield));
  return new Dirac_optimalDomainWall_EvenOdd(Dirac_node_,Kernel_.get(),Gfield,PauliVillars); 
}

/// Dirac_DWF4DfullSolv
DiracDWF4DfullFactory::DiracDWF4DfullFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"Kernel5d", MANDATORY);
  DiracFactory_.save(new DiracDWF5dFactory(node));
  XML::next_sibling(node, "SolverDWF", MANDATORY);
  SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
}

Dirac_optimalDomainWall_4D* DiracDWF4DfullFactory::
getDiracOperator4D(Field* const Gfield){

  DW5D_.save(DiracFactory_.get()->getDiracOperatorWL(Gfield));
  Fopr_.save(new Fopr_DdagD_Precondition(DW5D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));

  DW5D_PV_.save(DiracFactory_.get()->getDiracOperatorPV(Gfield));
  Fopr_PV_.save(new Fopr_DdagD_Precondition(DW5D_PV_.get()));
  Solver_PV_.save(SolverFactory_.get()->getSolver(Fopr_PV_.get()));

  return new Dirac_optimalDomainWall_4D_fullSolv(DW5D_.get(),DW5D_PV_.get(),
						 Solver_.get(),Solver_PV_.get());
}

/// Dirac_DWF4DeoSolv
DiracDWF4DeoFactory::DiracDWF4DeoFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"Kernel5d", MANDATORY);
  DiracFactory_.save(new DiracDWF5dFactory(node));
  DiracEOFactory_.save(new DiracDWF5dEvenOddFactory(node));
  XML::next_sibling(node,"SolverDWF", MANDATORY);
  SolverFactory_.save(SolverOperators::createSolverOperatorFactory(node));
}

Dirac_optimalDomainWall_4D* DiracDWF4DeoFactory::
getDiracOperator4D(Field* const Gfield){
  DW5D_.save(   DiracFactory_.get()->getDiracOperatorWL(Gfield));
  DW5D_EO_.save(DiracEOFactory_.get()->getDiracOperatorWL(Gfield));
  FoprEO_.save(new Fopr_DdagD(DW5D_EO_.get()));
  SolverEO_.save(SolverFactory_.get()->getSolver(FoprEO_.get()));
  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_.get(),
						  SolverEO_.get()));

  DW5DPV_.save(    DiracFactory_.get()->getDiracOperatorPV(Gfield));
  DW5D_EO_PV_.save(DiracEOFactory_.get()->getDiracOperatorPV(Gfield));
  FoprEO_PV_.save(new Fopr_DdagD(DW5D_EO_PV_.get()));
  SolverEO_PV_.save(SolverFactory_.get()->getSolver(FoprEO_PV_.get()));
  Inv_PV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_PV_.get(),
						     SolverEO_PV_.get()));
  
  return new Dirac_optimalDomainWall_4D_eoSolv(DW5D_.get(),DW5DPV_.get(),
					       Inv_.get(),Inv_PV_.get());
}

/// Dirac_staggered_EvenOdd
Dirac_staggered_EvenOdd* DiracStaggeredEvenOddFactory::
getDiracOperatorSTG(Dstagg::Dtype dt,Field* const Gfield){
  return new Dirac_staggered_EvenOdd(Dirac_node_,dt,Gfield); }

/// Dirac_staggered_EvenOdd_Adjoint
#if NC_==3
Dirac_staggered_EvenOdd_Adjoint* DiracStaggeredEvenOddAdjointFactory::
getDiracOperatorSTG(Dstagg::Dtype dt,Field* const Gfield){
  return new Dirac_staggered_EvenOdd_Adjoint(Dirac_node_,dt,Gfield); }
#endif
