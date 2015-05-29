/*! @file dirac_BFM_wrapper_factory.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 *
 * Time-stamp: <2015-05-15 16:42:07 cossu>
 */

#include "dirac_BFM_wrapper_factory.hpp"

/// DiracBFMoperator
DiracBFMoperatorFactory::DiracBFMoperatorFactory(const XML::node node)
 :Dirac_node_(node){
  XML::node current_node = Dirac_node_; 
  XML::descend(current_node, "Operator", MANDATORY);
  /// forces EO
  Dirac5D_EO_factory_.save(new DiracEvenOdd_DWF5dFactory(current_node));
}

Dirac_BFM_Wrapper* DiracBFMoperatorFactory::getDirac(const InputConfig& input) {
  return createDirac(input); }
Dirac_BFM_Wrapper* DiracBFMoperatorFactory::getDiracPV(const InputConfig& input) {
  return createDiracPV(input); } 

Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDirac(const InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDirac(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilsonEO_.get(), Regular5D);

}

Dirac_BFM_Wrapper* DiracBFMoperatorFactory::
createDiracPV(const InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDiracPV(Gfield));
  
  return new Dirac_BFM_Wrapper(Dirac_node_, 
			       &Gfield.gconf->data, 
			       DiracWilsonEO_.get(), PauliVillars5D);
}

