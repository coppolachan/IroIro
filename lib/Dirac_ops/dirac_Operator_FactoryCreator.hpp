/*! @file dirac_Operator_FactoryCreator.hpp 
 *  @brief Declaration of the global functions to generate Dirac_op factories
 * Time-stamp: <2014-07-04 11:23:40 noaki>
 */
#ifndef DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#define DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#include "dirac_Operator_Factory.hpp"

namespace Diracs {

  DiracFactory* createGeneralDiracFactory(const XML::node&);
  DiracWilsonLikeFactory* createDiracWilsonLikeFactory(const XML::node&);
  DiracWilsonLikeEvenOddFactory* createDiracWilsonLikeEvenOddFactory(const XML::node&);

  /// needed to abstract DWF4Deo and DWF4Dfull for qurak propagators
  DiracDWF4dFactory* createDiracDWF4dFactory(const XML::node&);

  /// needed to abstract DWF5D and DWF5dEvenOdd for HMC
  DiracDWF5dFactory* createDiracDWF5dFactory(const XML::node&);

  DiracWilsonLikeFactory* createGeneralDiracWilsonLikeFactory(const XML::node&);
  
  DiracStaggeredEvenOddLikeFactory* createDiracStaggeredEvenOddLikeFactory(const XML::node&);
}
#endif
