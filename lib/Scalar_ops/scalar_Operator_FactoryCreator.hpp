#ifndef SCALAR_OPERATOR_FACTORYCREATOR_INCLUDED
#define SCALAR_OPERATOR_FACTORYCREATOR_INCLUDED

#include "scalar_Operator_Factory.hpp"

namespace ScOps{
  
  ScalarOpFactory* createScalarOpFactory(const XML::node& node);
  void nullCheck(const XML::node&);
  void stopMsg(const XML::node&);
}
#endif
