/*! @file dirac_BFM_wrapper_factory_HDCG.cpp
 *  @brief Implementation of the FactoryCreator for HDCG operators
 *
 * Time-stamp: <2015-05-15 17:34:23 cossu>
 */

#include "dirac_BFM_DomainWall_4D_eo.hpp"
#include "dirac_BFM_wrapper_factory_HDCG.hpp"

///////////////////////// HDCG support
/// DiracBFM_HDCGoperator
DiracBFM_HDCGoperatorFactory::DiracBFM_HDCGoperatorFactory(const XML::node node)
 :Dirac_node_(node){
  XML::node current_node = Dirac_node_; 
  XML::descend(current_node, "Operator", MANDATORY);
  /// forces EO
  Dirac5D_EO_factory_.save(new DiracEvenOdd_DWF5dFactory(current_node));
}

Dirac_BFM_HDCG_Wrapper* DiracBFM_HDCGoperatorFactory::getDirac(const InputConfig& input) {
  return createDirac(input); }
Dirac_BFM_HDCG_Wrapper* DiracBFM_HDCGoperatorFactory::getDiracPV(const InputConfig& input) {
  return createDiracPV(input); } 

Dirac_BFM_HDCG_Wrapper* DiracBFM_HDCGoperatorFactory::
createDirac(const InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDirac(Gfield));
  
  return new Dirac_BFM_HDCG_Wrapper(Dirac_node_, 
				    &Gfield.gconf->data, 
				    DiracWilsonEO_.get(), Regular5D);
}

Dirac_BFM_HDCG_Wrapper* DiracBFM_HDCGoperatorFactory::
createDiracPV(const InputConfig& Gfield){
  DiracWilsonEO_.save(Dirac5D_EO_factory_.get()->getDiracPV(Gfield));
  
  return new Dirac_BFM_HDCG_Wrapper(Dirac_node_, 
				    &Gfield.gconf->data, 
				    DiracWilsonEO_.get(), PauliVillars5D);
}
