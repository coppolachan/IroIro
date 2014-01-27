/*! @file dirac_BFM_wrapper_factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2014-01-24 14:53:45 neo>
 */

#include "dirac_BFM_wrapper_factory.hpp"

/// DiracBFMoperator
Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDirac(InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDirac(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilsonEO_.get(), Regular5D);

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}

Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDiracPV(InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDiracPV(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilsonEO_.get(), PauliVillars5D);

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}


///////////////////// 4D Operator
DiracDWF4dBFMeoFactory::DiracDWF4dBFMeoFactory(XML::node node):Dirac_node_(node){
  node_BFM_ = Dirac_node_;
  XML::descend(node_BFM_,"Kernel5d",MANDATORY);
  DiracBFMFactory_.save(new DiracBFMoperatorFactory(node_BFM_));

  Solver_node_ = node;
  XML::descend(Solver_node_,"Solver_DWF-EO_BGQ",MANDATORY);
  SolverFactory_.save(new SolverCG_DWF_opt_Factory(Solver_node_));
}

Dirac_DomainWall_4D* DiracDWF4dBFMeoFactory::createDirac(InputConfig& input){
  BFM_Kernel_.save( DiracBFMFactory_.get()->getDirac(input));
  BFM_KernelPV_.save( DiracBFMFactory_.get()->getDiracPV(input));

  SolverEO_.save(SolverFactory_.get()->getSolver(BFM_Kernel_.get()));
  SolverEOpv_.save(SolverFactory_.get()->getSolver(BFM_KernelPV_.get()));

  BFM_Kernel_.get()->set_SolverParams(Solver_node_);
  BFM_KernelPV_.get()->set_SolverParams(Solver_node_);
  
  BFM_Kernel_.get()->initialize();
  BFM_KernelPV_.get()->initialize();


  Inv_.save(new EvenOddUtils::Inverter_WilsonLike(BFM_Kernel_.get()->getInternalEO(),SolverEO_.get()));
  InvPV_.save(new EvenOddUtils::Inverter_WilsonLike(BFM_KernelPV_.get()->getInternalEO(),SolverEOpv_.get()));




  XML::node current_node = node_BFM_;  
  XML::descend(current_node, "Operator", MANDATORY);  
  return new Dirac_DomainWall_4D_eoSolv(current_node,Inv_.get(),InvPV_.get());
}

