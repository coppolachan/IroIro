/*! @file dirac_BFM_wrapper_factory_HDCG.hpp 
 *  @brief Declaration of BFM Dirac operators factories
 *
 * Time-stamp: <2015-05-15 17:09:40 cossu>
 */
#ifndef DIRAC_BFM_HDCG_FACT_
#define DIRAC_BFM_HDCG_FACT_

#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
//#include "dirac_BFM_wrapper.hpp"
#include "dirac_BFM_HDCG.hpp"
#include "dirac_BFM_wrapper_factory.hpp"

/*! @brief Concrete class for creating Dirac DWF-5d operator with BFM with HDCG support */
class DiracBFM_HDCGoperatorFactory : public DiracDWF5dFactory {
  const XML::node Dirac_node_;

  // Internal operator must be EvenOdd
  RaiiFactoryObj<DiracEvenOdd_DWF5dFactory> Dirac5D_EO_factory_;
  RaiiFactoryObj<DiracWilsonLike_EvenOdd> DiracWilsonEO_;
  Dirac_BFM_HDCG_Wrapper* createDirac(const InputConfig&);
  Dirac_BFM_HDCG_Wrapper* createDiracPV(const InputConfig&);
public:
  DiracBFM_HDCGoperatorFactory(const XML::node node);
  Dirac_BFM_HDCG_Wrapper* getDirac(const InputConfig& );
  Dirac_BFM_HDCG_Wrapper* getDiracPV(const InputConfig& );
};


#endif
