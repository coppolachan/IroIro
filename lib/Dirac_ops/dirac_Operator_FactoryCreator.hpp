/*! @file dirac_Operator_FactoryCreator.hpp 
 *  @brief Declaration of the global functions to generate Dirac_op factories
 * Time-stamp: <2013-04-17 13:45:34 noaki>
 */
#ifndef DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#define DIRAC_OPERATOR_FACTORYCREATOR_INCLUDED
#include "dirac_Operator_Factory.hpp"

namespace DiracOperators {
  DiracOperatorFactory* 
  createGeneralDiracOperatorFactory(const XML::node);

  DiracWilsonLikeOperatorFactory* 
  createDiracWilsonLikeOperatorFactory(const XML::node);

  DiracDWF4dOperatorFactory* 
  createDiracDWF4dOperatorFactory(const XML::node);

  DiracDWF5dOperatorFactory* 
  createDiracDWF5dOperatorFactory(const XML::node);

  DiracWilsonLikeOperatorFactory* 
  createGeneralDiracWilsonLikeOperatorFactory(const XML::node);
  
  DiracStaggeredEvenOddLikeOperatorFactory* 
  createDiracStaggeredEvenOddLikeOperatorFactory(const XML::node);
}
#endif
