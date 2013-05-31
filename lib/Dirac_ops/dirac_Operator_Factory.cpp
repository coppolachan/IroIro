/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-05-27 17:05:39 noaki>
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
Dirac_Wilson* DiracWilsonFactory::getDiracWL(Field* const Gfield){
  return new Dirac_Wilson(Dirac_node_,Gfield);
}
/// Dirac_Wilson_EvenOdd
Dirac_Wilson_EvenOdd* DiracWilsonEvenOddFactory::
getDiracEO(Field* const Gfield){
  return new Dirac_Wilson_EvenOdd(Dirac_node_,Gfield);
}

/// Dirac_Wilson_Brillouin
Dirac_Wilson_Brillouin* DiracWilsonBrillouinFactory::
getDiracWL(Field* const Gfield){
  return new Dirac_Wilson_Brillouin(Dirac_node_,Gfield,type_);}

/// Dirac_Clover
Dirac_Clover* DiracCloverFactory::getDiracWL(Field* const Gfield){
  return new Dirac_Clover(Dirac_node_,Gfield); }

/// Dirac_Mobius
DiracMobiusFactory::DiracMobiusFactory(XML::node node)
  :Dirac_node_(node){
  XML::node dnode = node;
  XML::descend(dnode,"Dirac_denominator",MANDATORY);
  XML::descend(node,"Solver",MANDATORY);

  DiracFactory_.save(Diracs::createDiracWilsonLikeFactory(dnode));
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

Dirac_Mobius* DiracMobiusFactory::getDiracWL(Field* const Gfield){
  D_.save(DiracFactory_.get()->getDiracWL(Gfield));
  Fopr_.save(new Fopr_DdagD_Precondition(D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));
  return new Dirac_Mobius(Dirac_node_,D_.get(),Solver_.get());
}

/// Dirac_optimalDomainWall
DiracDomainWall5dFactory::DiracDomainWall5dFactory(XML::node node):Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}

Dirac_optimalDomainWall* DiracDomainWall5dFactory::getDiracWL(Field* const Gfield){
  return new Dirac_optimalDomainWall(Dirac_node_,
				     KernelFactory_.get()->getDiracWL(Gfield),
				     Gfield); 
}

Dirac_optimalDomainWall* DiracDomainWall5dFactory::getDiracPV(Field* const Gfield){
  return new Dirac_optimalDomainWall(Dirac_node_,
				     KernelFactory_.get()->getDiracWL(Gfield),
				     Gfield,PauliVillars); 
}

/// Dirac_optimalDomainWall_EvenOdd
DiracDomainWall5dEvenOddFactory::DiracDomainWall5dEvenOddFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeEvenOddFactory(node));
}

Dirac_optimalDomainWall_EvenOdd* 
DiracDomainWall5dEvenOddFactory::getDiracWL(Field* const Gfield){
  Kernel_.save(KernelFactory_.get()->getDiracEO(Gfield));
  return new Dirac_optimalDomainWall_EvenOdd(Dirac_node_,Kernel_.get(),Gfield); 
}

Dirac_optimalDomainWall_EvenOdd* 
DiracDomainWall5dEvenOddFactory::getDiracPV(Field* const Gfield){
  Kernel_.save(KernelFactory_.get()->getDiracEO(Gfield));
  return new Dirac_optimalDomainWall_EvenOdd(Dirac_node_,Kernel_.get(),
					     Gfield,PauliVillars); 
}

/// Dirac_DWF4DfullSolv
DiracDWF4DfullFactory::DiracDWF4DfullFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"Kernel5d", MANDATORY);
  DiracFactory_.save(new DiracDomainWall5dFactory(node));
  XML::next_sibling(node, "SolverDWF", MANDATORY);
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

Dirac_optimalDomainWall_4D* DiracDWF4DfullFactory::
getDirac4D(Field* const Gfield){

  DW5D_.save(DiracFactory_.get()->getDiracWL(Gfield));
  Fopr_.save(new Fopr_DdagD_Precondition(DW5D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));

  DW5D_PV_.save(DiracFactory_.get()->getDiracPV(Gfield));
  Fopr_PV_.save(new Fopr_DdagD_Precondition(DW5D_PV_.get()));
  Solver_PV_.save(SolverFactory_.get()->getSolver(Fopr_PV_.get()));

  return new Dirac_optimalDomainWall_4D_fullSolv(DW5D_.get(),DW5D_PV_.get(),
						 Solver_.get(),Solver_PV_.get());
}

/// Dirac_DWF4DeoSolv
DiracDWF4DeoFactory::DiracDWF4DeoFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"Kernel5d", MANDATORY);
  DiracEOFactory_.save(new DiracDomainWall5dEvenOddFactory(node));
  XML::next_sibling(node,"SolverDWF", MANDATORY);
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

Dirac_optimalDomainWall_4D* DiracDWF4DeoFactory::
getDirac4D(Field* const Gfield){
  DW5D_EO_.save(DiracEOFactory_.get()->getDiracWL(Gfield));
  FoprEO_.save(new Fopr_DdagD(DW5D_EO_.get()));
  SolverEO_.save(SolverFactory_.get()->getSolver(FoprEO_.get()));
  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_.get(),SolverEO_.get()));

  DW5D_EO_PV_.save(DiracEOFactory_.get()->getDiracPV(Gfield));
  FoprEO_PV_.save(new Fopr_DdagD(DW5D_EO_PV_.get()));
  SolverEO_PV_.save(SolverFactory_.get()->getSolver(FoprEO_PV_.get()));
  Inv_PV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5D_EO_PV_.get(),SolverEO_PV_.get()));
  
  return new Dirac_optimalDomainWall_4D_eoSolv(Dirac_node_,Inv_.get(),Inv_PV_.get());
}

/// Dirac_staggered_EvenOdd
Dirac_staggered_EvenOdd* DiracStaggeredEvenOddFactory::
getDiracSTG(Dstagg::Dtype dt,Field* const Gfield){
  return new Dirac_staggered_EvenOdd(Dirac_node_,dt,Gfield); }

/// Dirac_staggered_EvenOdd_Adjoint
#if NC_==3
Dirac_staggered_EvenOdd_Adjoint* DiracStaggeredEvenOddAdjointFactory::
getDiracSTG(Dstagg::Dtype dt,Field* const Gfield){
  return new Dirac_staggered_EvenOdd_Adjoint(Dirac_node_,dt,Gfield); }
#endif
