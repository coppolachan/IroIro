/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators

 * Time-stamp: <2014-07-04 17:48:58 noaki>
 */
#include "dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"
#include "Communicator/comm_io.hpp"
#include "eoUtils.hpp"
#include <string.h>

/// Dirac_Wilson
DiracWilsonLike* DiracWilsonFactory::createDirac(InputConfig& input){
  return new Dirac_Wilson(Dirac_node_,input.getGconf());}

/// Dirac_Wilson_Adjoint
DiracWilsonLike* DiracWilsonAdjointFactory::createDirac(InputConfig& input){
  return new Dirac_Wilson_Adjoint(Dirac_node_,input.getGconf());}

/// Dirac_Wilson_Brillouin
DiracWilsonLike* DiracWilsonBrillouinFactory::createDirac(InputConfig& input){
  return new Dirac_Wilson_Brillouin(Dirac_node_,input.getGconf(),type_);
}
/// Dirac_Wilson_Brillouin_OSS
DiracWilsonLike* DiracWilsonBrillouinOSSFactory::createDirac(InputConfig& input){
  return new Dirac_Wilson_Brillouin_OSS(Dirac_node_,input.getGconf(),type_);
}

/// Dirac_Wilson_FiniteDensity
DiracWilsonLike* DiracWilsonFiniteDensityFactory::createDirac(InputConfig& input){
  return new Dirac_Wilson_FiniteDensity(Dirac_node_,input.getGconf());
}

/// Dirac_Clover
DiracCloverFactory::DiracCloverFactory(XML::node node):Dirac_node_(node){
  XML::descend(node,"BaseWilson",MANDATORY);
  DiracFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}
DiracWilsonLike* DiracCloverFactory::createDirac(InputConfig& input){
  D_.save(DiracFactory_.get()->getDirac(input));
  return new Dirac_Clover(Dirac_node_,D_.get(),input.getGconf());}

/// Dirac_Mobius
DiracMobiusFactory::DiracMobiusFactory(XML::node node):Dirac_node_(node){
  XML::node dnode = node;
  XML::descend(dnode,"Dirac_denominator",MANDATORY);
  XML::descend(node,"Solver",MANDATORY);

  DiracFactory_.save(Diracs::createDiracWilsonLikeFactory(dnode));
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

DiracWilsonLike* DiracMobiusFactory::createDirac(InputConfig& input){
  D_.save(DiracFactory_.get()->getDirac(input));
  Fopr_.save(new Fopr_DdagD(D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));
  return new Dirac_Mobius(Dirac_node_,D_.get(),Solver_.get());
}

/// Dirac_DomainWall
DiracDomainWall5dFactory::DiracDomainWall5dFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}
Dirac_DomainWall* DiracDomainWall5dFactory::createDirac(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall(Dirac_node_,Kernel_.get()); 
}
Dirac_DomainWall* DiracDomainWall5dFactory::createDiracPV(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall(Dirac_node_,Kernel_.get(),PauliVillars); 
}

/// Dirac_DomainWall_Adjoint
DiracDomainWall5dAdjointFactory::DiracDomainWall5dAdjointFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}
Dirac_DomainWall_Adjoint* DiracDomainWall5dAdjointFactory::
createDirac(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_Adjoint(Dirac_node_,Kernel_.get()); 
}
Dirac_DomainWall_Adjoint* DiracDomainWall5dAdjointFactory::
createDiracPV(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_Adjoint(Dirac_node_,Kernel_.get(),PauliVillars); 
}

/// Dirac_Wilson_EvenOdd
DiracWilsonLike_EvenOdd* DiracWilsonEvenOddFactory::
createDirac(InputConfig& input){
  return new Dirac_Wilson_EvenOdd<Dirac_Wilson>(Dirac_node_,input.getGconf());}

/// Dirac_Wilson_Adjoint_EvenOdd
DiracWilsonLike_EvenOdd* DiracWilsonAdjointEvenOddFactory::
createDirac(InputConfig& input){
  return new Dirac_Wilson_EvenOdd<Dirac_Wilson_Adjoint>(Dirac_node_,input.getGconf());}

/// Dirac_DomainWall_EvenOdd
DiracEvenOdd_DWF5dFactory::DiracEvenOdd_DWF5dFactory(XML::node node)
:Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeEvenOddFactory(node));
}

Dirac_DomainWall_EvenOdd* DiracEvenOdd_DWF5dFactory::
createDirac(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_EvenOdd(Dirac_node_,Kernel_.get());
}
Dirac_DomainWall_EvenOdd* DiracEvenOdd_DWF5dFactory::
createDiracPV(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_EvenOdd(Dirac_node_,Kernel_.get(),
				      PauliVillars); 
}

/// Dirac_DomainWall_Adjoint_EvenOdd
DiracEvenOdd_DWF5dAdjointFactory::
DiracEvenOdd_DWF5dAdjointFactory(XML::node node):Dirac_node_(node){
  XML::descend(node,"BaseKernel", MANDATORY);
  KernelFactory_.save(Diracs::createDiracWilsonLikeEvenOddFactory(node));
}

Dirac_DomainWall_Adjoint_EvenOdd* DiracEvenOdd_DWF5dAdjointFactory::
createDirac(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_Adjoint_EvenOdd(Dirac_node_,Kernel_.get());
}
Dirac_DomainWall_Adjoint_EvenOdd* DiracEvenOdd_DWF5dAdjointFactory::
createDiracPV(InputConfig& input){
  if(!Kernel_.is_saved())
    Kernel_.save(KernelFactory_.get()->getDirac(input));
  return new Dirac_DomainWall_Adjoint_EvenOdd(Dirac_node_,Kernel_.get(),
					      PauliVillars); 
}

/// Dirac_DWF4DfullSolv
DiracDWF4DfullFactory::DiracDWF4DfullFactory(XML::node node,DW5dPrecond prec)
:Dirac_node_(node),prec_(prec){
  XML::descend(node,"Kernel5d", MANDATORY);
  DiracFactory_.save(new DiracDomainWall5dFactory(node));
  XML::next_sibling(node,"SolverDWF", MANDATORY);
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

Dirac_DomainWall_4D* DiracDWF4DfullFactory::createDirac(InputConfig& input){
  DW5D_.save(DiracFactory_.get()->getDirac(input));
  Fopr_.save(new Fopr_DdagD(DW5D_.get()));
  Solver_.save(SolverFactory_.get()->getSolver(Fopr_.get()));

  DW5dPV_.save(new Dirac_DomainWall(*(DW5D_.get()),PauliVillars));
  FoprPV_.save(new Fopr_DdagD(DW5dPV_.get()));
  SolverPV_.save(SolverFactory_.get()->getSolver(FoprPV_.get()));

  return new Dirac_DomainWall_4D_fullSolv(DW5D_.get(),DW5dPV_.get(),
						 Solver_.get(),SolverPV_.get(),prec_);
}

/// Dirac_DWF4DeoSolv
DiracDWF4DeoFactory::DiracDWF4DeoFactory(XML::node node):Dirac_node_(node){
  XML::descend(Dirac_node_,"Kernel5d",MANDATORY);

  XML::node node_DWF5d = Dirac_node_;
  DiracEOFactory_.save(new DiracEvenOdd_DWF5dFactory(node_DWF5d));
  XML::descend(node,"SolverDWF",MANDATORY);
  SolverFactory_.save(Solvers::createSolverFactory(node));
}

Dirac_DomainWall_4D* DiracDWF4DeoFactory::createDirac(InputConfig& input){
  DW5dEO_.save(DiracEOFactory_.get()->getDirac(input));
  FoprEO_.save(new Fopr_DdagD(DW5dEO_.get()));
  SolverEO_.save(SolverFactory_.get()->getSolver(FoprEO_.get()));
  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5dEO_.get(),SolverEO_.get()));

  //DW5dEOpv_.save(new Dirac_DomainWall_EvenOdd(*(DW5dEO_.get()),PauliVillars));
  DW5dEOpv_.save(DiracEOFactory_.get()->getDiracPV(input));
  FoprEOpv_.save(new Fopr_DdagD(DW5dEOpv_.get()));
  SolverEOpv_.save(SolverFactory_.get()->getSolver(FoprEOpv_.get()));
  InvPV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5dEOpv_.get(),SolverEOpv_.get()));

  return new Dirac_DomainWall_4D_eoSolv(Dirac_node_,Inv_.get(),InvPV_.get());
}

#ifdef IBM_BGQ_WILSON
/// Dirac_DWF4dBGQeoSolv
DiracDWF4dBGQeoFactory::DiracDWF4dBGQeoFactory(XML::node node):Dirac_node_(node){
  XML::descend(Dirac_node_,"Kernel5d",MANDATORY);
  XML::node node_DWF5d = Dirac_node_;
  DiracEOFactory_.save(new DiracEvenOdd_DWF5dFactory(node_DWF5d));

  XML::descend(node,"SolverDWF",MANDATORY);
  SolverFactory_.save(new SolverCG_DWF_opt_Factory(node));
}

Dirac_DomainWall_4D* DiracDWF4dBGQeoFactory::createDirac(InputConfig& input){
  DW5dEO_.save(DiracEOFactory_.get()->getDirac(input));
  DW5dEOpv_.save(new Dirac_DomainWall_EvenOdd(*(DW5dEO_.get()),PauliVillars));

  SolverEO_.save(SolverFactory_.get()->getSolver(DW5dEO_.get()));
  SolverEOpv_.save(SolverFactory_.get()->getSolver(DW5dEOpv_.get()));

  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(DW5dEO_.get(),SolverEO_.get()));
  InvPV_.save(new EvenOddUtils::Inverter_WilsonLike(DW5dEOpv_.get(),SolverEOpv_.get()));

  return new Dirac_DomainWall_4D_eoSolv(Dirac_node_,Inv_.get(),InvPV_.get());
}
#endif

/// Dirac_DWoverlap
DiracDWoverlapFactory::DiracDWoverlapFactory(XML::node node):Dirac_node_(node){
  XML::descend(node,"Kernel4d",MANDATORY);
  DW4dFactory_.save(Diracs::createDiracDWF4dFactory(node));
}

DiracWilsonLike* DiracDWoverlapFactory::createDirac(InputConfig& input){
  DW4d_.save(DW4dFactory_.get()->getDirac(input));
  return new Dirac_DWoverlap(DW4d_.get(),input.emode);
}

/// Dirac_LowModeDeflation with exact eigenmodes
DiracDeflationExactFactory::DiracDeflationExactFactory(XML::node node)
 :Dirac_node_(node){
  XML::descend(node,"KernelCore",MANDATORY);
  DwFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}

Dirac_LowModeDeflation* DiracDeflationExactFactory::createDirac(InputConfig& input){
  Dw_.save(DwFactory_.get()->getDirac(input));
  return new Dirac_LowModeDeflation_ExactEigen(Dw_.get(),input.emode);
}

/// Dirac_LowModeDeflation with approximated subspace
DiracDeflationApproxFactory::DiracDeflationApproxFactory(XML::node node)
 :Dirac_node_(node){
  XML::descend(node,"KernelCore",MANDATORY);
  DwFactory_.save(Diracs::createDiracWilsonLikeFactory(node));
}

Dirac_LowModeDeflation* DiracDeflationApproxFactory::createDirac(InputConfig& input){
  Dw_.save(DwFactory_.get()->getDirac(input));
  return new Dirac_LowModeDeflation_Approx(Dw_.get(),input.emode);
}

/// Dirac_staggered_EvenOdd
Dirac_staggered_EvenOdd* DiracStaggeredEvenOddFactory::createDirac(InputConfig& input){
  return new Dirac_staggered_EvenOdd(Dirac_node_,input.getGconf(),Dstagg::DdagDee); 
}
Dirac_staggered_EvenOdd* DiracStaggeredEvenOddFactory::createDoo(InputConfig& input){
  return new Dirac_staggered_EvenOdd(Dirac_node_,input.getGconf(),Dstagg::DdagDoo); 
}

/// Dirac_staggered_EvenOdd_Adjoint
#if NC_==3
Dirac_staggered_EvenOdd_Adjoint* DiracStaggeredEvenOddAdjointFactory::
createDirac(InputConfig& input){
  return new Dirac_staggered_EvenOdd_Adjoint(Dirac_node_,input.getGconf(),
					     Dstagg::DdagDee); 
}
Dirac_staggered_EvenOdd_Adjoint* DiracStaggeredEvenOddAdjointFactory::
createDoo(InputConfig& input){
  return new Dirac_staggered_EvenOdd_Adjoint(Dirac_node_,input.getGconf(),
					     Dstagg::DdagDoo); 
}
#endif
