#include "scalar_Operator_FactoryCreator.hpp"
#include "PugiXML/xmlUtilities.hpp"

namespace ScOps{
  ScalarOpFactory* createScalarOpFactory(const XML::node& node){
    XML::nullCheck(node,"ScalarOp");
    const char* Scalar_name = node.attribute("name").value();
    
    if(!strcmp(Scalar_name,"Laplacian"))   return new LaplacianFactory(node);
    if(!strcmp(Scalar_name,"LaplacianXt")) return new Laplacian4DSfactory();
    if(!strcmp(Scalar_name,"LaplacianXtXd")) return new Laplacian4DFfactory();
    XML::stopMsg(node,"ScalarOp"); 
  }
}
