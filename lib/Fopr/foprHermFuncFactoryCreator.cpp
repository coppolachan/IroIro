#include "foprHermFuncFactoryCreator.hpp"
#include "PugiXML/xmlUtilities.hpp"

namespace FuncHermite{
  
  FoprHermFactory* createHermOpFuncFactory(const XML::node& node,
					   std::string str){
    XML::nullCheck(node,str.c_str());
    const char* Fname = node.attribute("name").value();

    if(!strcmp(Fname,"NULL"))       return new FoprNULLfuncFactory();
    if(!strcmp(Fname,"Chebyshev"))  return new FoprChebyshevFuncFactory(node);
    if(!strcmp(Fname,"Linear"))     return new FoprLinearFuncFactory(node);
    if(!strcmp(Fname,"QuadLinear")) return new FoprQuadLinearFuncFactory(node);
    if(!strcmp(Fname,"Exponential"))return new FoprExpFuncFactory(node);
    if(!strcmp(Fname, "LimitExp"))  return new FoprLexpFuncFactory(node);
    
    XML::stopMsg(node,str.c_str());
  }
}// end of namespace
