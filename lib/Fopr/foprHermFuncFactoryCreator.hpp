#ifndef FOPRHERMFUNCFACTORYCREATOR_INCLUDED
#define FOPRHERMFUNCFACTORYCREATOR_INCLUDED
#include "foprHermFuncFactory.hpp"
#include <string>

namespace FuncHermite{
  FoprHermFactory* createHermOpFuncFactory(const XML::node&,
					   std::string str="HermiteOp");
}
#endif
