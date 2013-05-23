/*! @file dirac_Operator_FactoryCreator.hpp 
 *  @brief Declaration of the global functions to generate Dirac_op factories
 * Time-stamp: <2013-05-23 11:14:05 noaki>
 */
#ifndef DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#define DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#include "dirac_Operator_Factory.hpp"

namespace Diracs {

  DiracFactory* createGeneralDiracFactory(const XML::node&);
  DiracWilsonLikeFactory* createDiracWilsonLikeFactory(const XML::node&);

  DiracWilsonLikeEvenOddFactory* 
  createDiracWilsonLikeEvenOddFactory(const XML::node&);

  DiracDWF4dFactory* createDiracDWF4dFactory(const XML::node&);

  DiracDWF5dFactory* createDiracDWF5dFactory(const XML::node&);

  DiracWilsonLikeFactory* 
  createGeneralDiracWilsonLikeFactory(const XML::node&);
  
  DiracStaggeredEvenOddLikeFactory* 
  createDiracStaggeredEvenOddLikeFactory(const XML::node&);

  void nullCheck(const XML::node&);
  void stopMsg(const XML::node&);
}
#endif
