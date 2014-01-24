#include "chebyshevAccelFunc.hpp"
#include "pugi_interface.h"
#include <string.h>

ChebyshevAccelFunc::ChebyshevAccelFunc(XML::node node){

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
    mapFact_.reset(new FoprLinearFunc(slope,itcpt)); 
    
  }else if(!strcmp(kn_name,"QuadLinear")){
    slope = sgn*2.0/(exu*exu-exl*exl);
    itcpt = -sgn*(exu*exu+exl*exl)/(exu*exu-exl*exl);
    mapFact_.reset(new FoprQuadLinearFunc(slope,itcpt)); 
    
  }else{
    CCIO::cout<<kn_name<<" is not compatible with current implementation.\n";
    abort();
  }
    
  //// chebyshev function ////
  XML::descend(node,"ChebyshevFunc");
  cbfunc_.reset(new FoprChebyshevFunc(node));
}							  
  
Fopr_Chebyshev* ChebyshevAccelFunc::getFoprHerm(const Fopr_Herm* foprH)const{
  mapOp_.reset(mapFact_.get()->getFoprHerm(foprH));
  return cbfunc_.get()->getFoprHerm(mapOp_.get());
}



