/*! @file dirac_Operator_Factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-07-12 17:34:18 cossu>
 */

#include "dirac_BFM_wrapper_factory.hpp"

/// DiracBFMoperator
Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
getDiracWL(Field* const Gfield){
  DiracWilson_.save(DiracWL_factory_.get()->getDiracWL(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       Gfield, 
			       DiracWilson_.get());

  //if (XML::descend(Dirac_node_, "Solver"))
  //  BFMop_.get()->set_SolverParams(Dirac_node_));
  
}
