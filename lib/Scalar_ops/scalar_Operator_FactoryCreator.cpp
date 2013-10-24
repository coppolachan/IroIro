#include "scalar_Operator_FactoryCreator.hpp"

namespace ScOps{
  ScalarOpFactory* createScalarOpFactory(const XML::node& node){
    nullCheck(node);
    const char* Scalar_name = node.attribute("name").value();
    
    if(!strcmp(Scalar_name,"Laplacian"))
      return new LaplacianFactory(node);
    stopMsg(node);
  }

  void nullCheck(const XML::node& node){
    if(node==NULL){
      std::cout<<"Mandatory node is missing in input file(Scalar Obj)\n";
      abort();
    }else{
      const char* Scalar_name = node.attribute("name").value();
      if(!strcmp(Scalar_name,"")){
	std::cerr<<"No name provided for Scalar Operator. Request by <"
		 << node.name()<<">\n";
	abort();
      }
    }
  }

  void stopMsg(const XML::node& node){
    const char* Scalar_name = node.attribute("name").value();
    std::cerr<<"No Scalar Operator available with name ["
             << Scalar_name << "]. Request by <"<< node.name()<<">\n";
    abort();
  }
}
