#include "foprHermFactoryCreator.hpp"
#include "foprHermFuncFactoryCreator.hpp"
#include "PugiXML/xmlUtilities.hpp"

namespace HermiteOp{

  FoprHermFactory* createFoprHermFactory(XML::node node){
    XML::nullCheck(node,"HermiteianOperator");

    const char* hf_name = node.attribute("name").value();
    
    //// Hermitian Op. from scalar Op.
    if(!strcmp(hf_name,"Scalar")) return new FoprHermFactory_Scalar(node);

    //// Hermitian Op. as a functional of a Hermitian Op.
    if(!strcmp(hf_name,"Functional")) 
      return FuncHermite::createHermOpFuncFactory(node);

    //// Hermitian Op. from WilsonLikeDirac Op.
    if(!strcmp(hf_name,"HermD")) return new FoprHermFactory_HD(node);
    if(!strcmp(hf_name,"g5D"  )) return new FoprHermFactory_H(node);
    if(!strcmp(hf_name,"R5g5D")) return new FoprHermFactory_H5d(node);
    if(!strcmp(hf_name,"DdagD")) return new FoprHermFactory_DdagD(node);
    if(!strcmp(hf_name,"DDdag")) return new FoprHermFactory_DDdag(node);

    //// error msg.
    XML::stopMsg(node,"HermitianOperator");
    abort();
  }

}
