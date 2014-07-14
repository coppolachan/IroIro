#ifndef XMLUTILITIES_INCLUDED
#define XMLUTILITIES_INCLUDED
#include "pugi_interface.h"

namespace XML{
  void nullCheck(const XML::node& node,std::string ObjName);
  void stopMsg(const XML::node& node,std::string ObjName);
}

#endif
