/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-07-17 14:20:58 cossu>
 */

#include "dirac_BFM_wrapper_factory.hpp"

/// DiracBFMoperator
Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDirac(InputConfig& Gfield){
  DiracWilson_.save(Dirac5D_factory_.get()->getDirac(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilson_.get(), Regular5D);

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}

Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDiracPV(InputConfig& Gfield){
  DiracWilson_.save(Dirac5D_factory_.get()->getDiracPV(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilson_.get(), PauliVillars5D);

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}

