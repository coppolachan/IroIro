/*! @file dirac_Operator_FactoryCreator.cpp
 *  @brief Implementation of the FactoryCreator for Dirac operators
 * Time-stamp: <2014-07-04 11:29:56 noaki>
 */
#include "dirac_Operator_FactoryCreator.hpp"
#include "PugiXML/xmlUtilities.hpp"
#ifdef HAVE_LIBBFM
#include "./BFM_Wrapper/dirac_BFM_wrapper_factory.hpp"
#endif

namespace Diracs {

  DiracWilsonLikeFactory* 
  createDiracWilsonLikeFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_name");
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
    if(!strcmp(Dirac_name, "DiracWilson_Brillouin_OSS")) 
      return new DiracWilsonBrillouinOSSFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Brillouin_Imp_OSS")) 
      return new DiracWilsonBrillouinOSSFactory(node,Improved);
    if(!strcmp(Dirac_name, "DiracWilson_FiniteDensity")) 
      return new DiracWilsonFiniteDensityFactory(node);
    if(!strcmp(Dirac_name, "DiracClover"))  
      return new DiracCloverFactory(node);
    if(!strcmp(Dirac_name, "DiracMobius"))  
      return new DiracMobiusFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dAdjoint"))  
      return new DiracDomainWall5dAdjointFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
      return new DiracEvenOdd_DWF5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dAdjointEvenOdd"))  
      return new DiracEvenOdd_DWF5dAdjointFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
      return new DiracDWF4DfullFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_LUprecond"))  
      return new DiracDWF4DfullFactory(node,LUprecond);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_eo"))  
      return new DiracDWF4DeoFactory(node);
#ifdef IBM_BGQ_WILSON
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BGQeo"))  
      return new DiracDWF4dBGQeoFactory(node);
#ifdef HAVE_LIBBFM
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BFMeo"))  
      return new DiracDWF4dBFMeoFactory(node);
#endif
#endif
    if(!strcmp(Dirac_name, "DiracDWoverlap"))  
      return new DiracDWoverlapFactory(node);
    if(!strcmp(Dirac_name, "DiracDeflation_ExactEigen"))  
      return new DiracDeflationExactFactory(node);
    if(!strcmp(Dirac_name, "DiracDeflation_ApproxSubspace"))  
      return new DiracDeflationApproxFactory(node);
    XML::stopMsg(node,"Dirac_name");
  }

  DiracWilsonLikeEvenOddFactory* 
  createDiracWilsonLikeEvenOddFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_name");
    const char* Dirac_name = node.attribute("name").value();

    if(!strcmp(Dirac_name, "DiracWilson"))
      return new DiracWilsonEvenOddFactory(node);
    if(!strcmp(Dirac_name, "DiracWilson_Adjoint"))
      return new DiracWilsonAdjointEvenOddFactory(node);

    XML::stopMsg(node,"Dirac_name");
  }

  DiracDWF4dFactory* createDiracDWF4dFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_name");
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
#ifdef HAVE_LIBBFM
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BFMeo"))  
      return new DiracDWF4dBFMeoFactory(node);
#endif
#endif
    XML::stopMsg(node,"Dirac_name");
  }

  // this one forces to have a PauliVillars operator
  // useful in the case of the DomainWall5d action
  DiracDWF5dFactory* createDiracDWF5dFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_name");
    const char* Dirac_name = node.attribute("name").value();

    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dAdjoint"))  
      return new DiracDomainWall5dAdjointFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dEvenOdd"))  
      return new DiracEvenOdd_DWF5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5dAdjointEvenOdd"))  
      return new DiracEvenOdd_DWF5dAdjointFactory(node);
    XML::stopMsg(node,"Dirac_name");      
  }

  // This is for the Dirac-op creation beyond the WilsonLike framework
  DiracFactory* createGeneralDiracFactory(const XML::node& node){
    // temporal hack
    return createGeneralDiracWilsonLikeFactory(node);
  }

  DiracWilsonLikeFactory* 
  createGeneralDiracWilsonLikeFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_name");
    const char* Dirac_name = node.attribute("name").value();
    
    if(!strcmp(Dirac_name, "DiracWilson"))  
      return new DiracWilsonFactory(node);
    if(!strcmp(Dirac_name, "DiracWilonFiniteDensity"))  
      return new DiracWilsonFiniteDensityFactory(node);
    if(!strcmp(Dirac_name, "DiracClover"))  
      return new DiracCloverFactory(node);
    if(!strcmp(Dirac_name, "DiracMobius"))  
      return new DiracMobiusFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall5d"))  
      return new DiracDomainWall5dFactory(node);
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d"))  
      return new DiracDWF4DfullFactory(node);
#ifdef IBM_BGQ_WILSON
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BGQeo"))  
      return new DiracDWF4dBGQeoFactory(node);
#ifdef HAVE_LIBBFM
    if(!strcmp(Dirac_name, "DiracOptimalDomainWall4d_BFMeo"))  
      return new DiracDWF4dBFMeoFactory(node);
#endif
#endif
    XML::stopMsg(node,"Dirac_name");      
  }

  DiracStaggeredEvenOddLikeFactory* 
  createDiracStaggeredEvenOddLikeFactory(const XML::node& node){
    XML::nullCheck(node,"Dirac_node");
    const char* Dirac_name = node.attribute("name").value();
    
    if(!strcmp(Dirac_name, "DiracStaggeredEvenOdd"))  
      return new DiracStaggeredEvenOddFactory(node);
    
#if NC_==3 //Only for NC=3
    if (!strcmp(Dirac_name, "DiracStaggeredEvenOddAdjoint"))  
      return new DiracStaggeredEvenOddAdjointFactory(node);
#endif
    XML::stopMsg(node,"Dirac_name");      
  }

} //end of namespace
