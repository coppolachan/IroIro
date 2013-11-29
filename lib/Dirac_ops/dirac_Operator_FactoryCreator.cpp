/*! @file dirac_Operator_FactoryCreator.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2013-11-29 16:36:17 noaki>
 */
#include "dirac_Operator_FactoryCreator.hpp"

namespace Diracs {

  DiracWilsonLikeFactory* 
  createDiracWilsonLikeFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();
    
    if(!strcmp(Dirac_name, "DiracWilson"))  
      return new DiracWilsonFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Adjoint"))  
      return new DiracWilsonAdjointFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_EvenOdd")) 
      return new DiracWilsonEvenOddFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Adjoint_EvenOdd")) 
      return new DiracWilsonAdjointEvenOddFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Brillouin")) 
      return new DiracWilsonBrillouinFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Brillouin_Imp")) 
      return new DiracWilsonBrillouinFactory(node,Improved);
    if(!strcmp(Dirac_name, "DiracClover"))  
      return new DiracCloverFactory(node);
    if(!strcmp(Dirac_name, "DiracMobius"))  
      return new DiracMobiusFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
      return new DiracEvenOdd_DWF5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
      return new DiracDWF4DfullFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_LUprecond"))  
      return new DiracDWF4DfullFactory(node,LUprecond);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_eo"))  
      return new DiracDWF4DeoFactory(node);
#ifdef IBM_BGQ_WILSON
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BGQeo"))  
      return new DiracDWF4dBGQeoFactory(node);
#endif
    if(!strcmp(Dirac_name, "DiracDWoverlap"))  
      return new DiracDWoverlapFactory(node);
    if(!strcmp(Dirac_name, "DiracDeflation_ExactEigen"))  
      return new DiracDeflationExactFactory(node);
    if(!strcmp(Dirac_name, "DiracDeflation_ApproxSubspace"))  
      return new DiracDeflationApproxFactory(node);
    stopMsg(node);
  }

  DiracWilsonLikeEvenOddFactory* 
  createDiracWilsonLikeEvenOddFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();

    if(!strcmp(Dirac_name, "DiracWilson"))
      return new DiracWilsonEvenOddFactory(node);
    if(!strcmp(Dirac_name, "DiracWilsonAdjoint"))
      return new DiracWilsonAdjointEvenOddFactory(node);

    stopMsg(node);
  }

  DiracDWF4dFactory* createDiracDWF4dFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();

    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
      return new DiracDWF4DfullFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_LUprecond"))  
      return new DiracDWF4DfullFactory(node,LUprecond);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_eo"))  
      return new DiracDWF4DeoFactory(node);
#ifdef IBM_BGQ_WILSON
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BGQeo"))  
      return new DiracDWF4dBGQeoFactory(node);
#endif
    stopMsg(node);
  }

  // this one forces to have a PauliVillars operator
  // useful in the case of the DomainWall5d action
  DiracDWF5dFactory* createDiracDWF5dFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();

    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
      return new DiracEvenOdd_DWF5dFactory(node);
    stopMsg(node);      
  }

  /*
  // This is for the Dirac-op creation beyond the WilsonLike framework
  DiracFactory* createGeneralDiracFactory(const XML::node& node){
    // temporal hack
    return createGeneralDiracWilsonLikeFactory(node);
  }

  DiracWilsonLikeFactory* 
  createGeneralDiracWilsonLikeFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();
    
    if(!strcmp(Dirac_name, "DiracWilson"))  
      return new DiracWilsonFactory(node);
    if(!strcmp(Dirac_name, "DiracClover"))  
      return new DiracCloverFactory(node);
    if(!strcmp(Dirac_name, "DiracMobius"))  
      return new DiracMobiusFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
      return new DiracDWF4DfullFactory(node);
    stopMsg(node);      
  }
  */

  DiracStaggeredEvenOddLikeFactory* 
  createDiracStaggeredEvenOddLikeFactory(const XML::node& node){
    nullCheck(node);
    const char* Dirac_name = node.attribute("name").value();
    
    if(!strcmp(Dirac_name, "DiracStaggeredEvenOdd"))  
      return new DiracStaggeredEvenOddFactory(node);
    
#if NC_==3 //Only for NC=3
    if (!strcmp(Dirac_name, "DiracStaggeredEvenOddAdjoint"))  
      return new DiracStaggeredEvenOddAdjointFactory(node);
#endif
    stopMsg(node);      
  }

  //////////////// utility functions ////////////////////
  void nullCheck(const XML::node& node){
    if(node==NULL){
      std::cout<<"Mandatory node is missing in input file(Dirac obj)\n";
      abort();
    }else{
      const char* Dirac_name = node.attribute("name").value();
      if(!strcmp(Dirac_name,"")){
	std::cerr<<"No name provided for Dirac Operator. Request by <"
		 << node.name()<<">\n";
	abort();
      }
    }
  }

  void stopMsg(const XML::node& node){
    const char* Dirac_name = node.attribute("name").value();
    std::cerr<<"No Dirac Operator available with name ["
             << Dirac_name << "]. Request by <"<< node.name()<<">\n";
    abort();
  }

} //end of namespace
