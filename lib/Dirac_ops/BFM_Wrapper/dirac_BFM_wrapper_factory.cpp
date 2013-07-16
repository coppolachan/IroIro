/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-07-12 18:37:28 cossu>
 */

#include "dirac_BFM_wrapper_factory.hpp"

/// DiracBFMoperator
Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDirac(InputConfig& Gfield){
  DiracWilson_.save(DiracWL_factory_.get()->getDirac(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilson_.get());

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}

