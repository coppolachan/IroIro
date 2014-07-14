#include "chebyshevAccelFuncFactory.hpp"
#include "pugi_interface.h"
#include "PugiXML/xmlUtilities.hpp"
#include <string.h>

ChebyshevAccelFuncFactory::
ChebyshevAccelFuncFactory(XML::node node){

  double sgn = 1.0;
  const char* st_name = node.attribute("sorting").value();
  if(!strcmp(st_name,"Lowest" )) sgn = -1.0;
  
  double thrs;
  XML::read(node,"threshold",thrs);
  XML::descend(node,"Acceleration");
  
  //// mapping function ////
  XML::node mp_node = node;
  XML::descend(mp_node,"Mapping");
  const char* kn_name = mp_node.attribute("function").value();

  double exl,exu;
  XML::read(mp_node,"ex_lower",exl,MANDATORY);
  if( sgn < 0.0 && exl < thrs){
    CCIO::cout<<"it must be threshold < ex_lower for the lowest sorting.\n";
    abort();
  }
  XML::read(mp_node,"ex_upper",exu,MANDATORY);
  if( sgn > 0.0 && exu > thrs){
    CCIO::cout<<"it must be threshold > ex_upper for the highest sorting.\n";
    abort();
  }

  double slope,itcpt;
  if(      !strcmp(kn_name,"Linear")){
    slope = sgn*2.0/(exu-exl);
    itcpt = -sgn*(exu+exl)/(exu-exl);
    mapFact_.reset(new LinearFuncFactory(slope,itcpt)); 

  }else if(!strcmp(kn_name,"QuadLinear")){
    slope = sgn*2.0/(exu*exu-exl*exl);
    itcpt = -sgn*(exu*exu+exl*exl)/(exu*exu-exl*exl);
    mapFact_.reset(new QuadLinearFuncFactory(slope,itcpt)); 
    
  }else{ XML::stopMsg(node,kn_name);}
    
  XML::descend(node,"ChebyshevFunc");
  cbNode_= node;
}							  
  
Fopr_Chebyshev* ChebyshevAccelFuncFactory::getOp(const Fopr_Herm* op)const{
  mapOp_.reset(mapFact_->getOp(op));
  return new Fopr_Chebyshev(cbNode_,mapOp_.get());
}


