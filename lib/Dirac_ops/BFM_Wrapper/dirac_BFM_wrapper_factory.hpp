/*! @file dirac_BFM_wrapper_factory.hpp 
 *  @brief Declaration of BFM Dirac operators factories
 Time-stamp: <2013-07-16 10:01:54 cossu>
 */
#ifndef DIRAC_BFM_FACT_
#define DIRAC_BFM_FACT_

#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "dirac_BFM_wrapper.hpp"

/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracBFMoperatorFactory : public DiracWilsonLikeFactory {
  const XML::node Dirac_node_;
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFMop_;

  RaiiFactoryObj<DiracWilsonLikeFactory> DiracWL_factory_;
  RaiiFactoryObj<DiracWilsonLike> DiracWilson_;
  Dirac_BFM_Wrapper* createDirac(InputConfig&);
public:
  DiracBFMoperatorFactory(const XML::node node):Dirac_node_(node){
    CCIO::cout << "Creating the factory\n";
    XML::node current_node = Dirac_node_; 
    XML::descend(current_node, "Operator", MANDATORY);
    DiracWL_factory_.save(Diracs::createDiracWilsonLikeFactory(current_node));
  
  }

  Dirac_BFM_Wrapper* getDirac(InputConfig& input) {
    return createDirac(input);
  }
  
};


#endif
