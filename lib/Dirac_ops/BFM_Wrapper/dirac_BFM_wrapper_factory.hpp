/*! @file dirac_BFM_wrapper_factory.hpp 
 *  @brief Declaration of BFM Dirac operators factories
 Time-stamp: <2013-07-17 14:20:44 cossu>
 */
#ifndef DIRAC_BFM_FACT_
#define DIRAC_BFM_FACT_

#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "dirac_BFM_wrapper.hpp"

/*! @brief Concrete class for creating Dirac Optimal DWF-5d operators */
class DiracBFMoperatorFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;
  RaiiFactoryObj<Dirac_BFM_Wrapper> BFMop_;

  RaiiFactoryObj<DiracDWF5dFactory> Dirac5D_factory_;
  RaiiFactoryObj<DiracWilsonLike> DiracWilson_;
  Dirac_BFM_Wrapper* createDirac(InputConfig&);
  Dirac_BFM_Wrapper* createDiracPV(InputConfig&);
public:
  DiracBFMoperatorFactory(const XML::node node):Dirac_node_(node){
    XML::node current_node = Dirac_node_; 
    XML::descend(current_node, "Operator", MANDATORY);
    Dirac5D_factory_.save(Diracs::createDiracDWF5dFactory(current_node));
  
  }

  Dirac_BFM_Wrapper* getDirac(InputConfig& input) {
    return createDirac(input);
  }
  Dirac_BFM_Wrapper* getDiracPV(InputConfig& input) {
    return createDiracPV(input);
  }
  
};


#endif
