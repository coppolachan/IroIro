/*!@file fopr_ChebyshevDdagDLin_Factory.hpp
 * @brief factory of fopr_Chebyshev whose kernel is DdagDLin
 This is especially used for the chebyshev acceleration in the eigenmodes calc.
 */
#ifndef FOPRHERMFACTORY_CHEBYSHEVDDAGDLIN_INCLUDED
#define FOPRHERMFACTORY_CHEBYSHEVDDAGDLIN_INCLUDED

#include "include/foprHermFactory.hpp"
#include "include/fopr_Chebyshev.h"

class FoprHermFactory_ChebyshevDdagDLin: public FoprHermFactory{
  int N_;
  FoprHermFactory_DdagDLinear ddlinFact_;
  std::auto_ptr<Fopr_Herm> kernel_;

public:
  FoprHermFactory_ChebyshevDdagDLin(const XML::node& node):kernel_(NULL){

    double sgn = 1.0;
    const char* st_name = node.attribute("sorting").value();
    if(!strcmp(st_name,"Lowest" )) sgn = -1.0;

    double thrs;
    XML::read(node,"threshold",thrs);

    XML::node ac_node = node;
    XML::descend(ac_node,"Acceleration");

    const char* kn_name = ac_node.attribute("kernel_op").value();
    
    if(!strcmp(kn_name,"DdagDLinear")){
      double exl,exu;
    
      XML::read(ac_node,"ex_lower",exl,MANDATORY);
      XML::read(ac_node,"ex_upper",exu,MANDATORY);
      XML::read(ac_node,"Npoly",N_,MANDATORY);
    
      if(!strcmp(st_name,"Lowest" ) && exl < thrs){
	CCIO::cout<<"it must be threshold < ex_lower for the lowest sorting.\n";
	abort();
      }
      if(!strcmp(st_name,"Highest" ) && exu > thrs){
	CCIO::cout<<"it must be threshold > ex_upper for the highest sorting.\n";
	abort();
      }

      double slope = sgn*2.0/(exu*exu-exl*exl);
      double itcpt = -sgn*(exu*exu+exl*exl)/(exu*exu-exl*exl);
    
      ddlinFact_= FoprHermFactory_DdagDLinear(slope,itcpt);

    }else{
      CCIO::cout<<kn_name<<" is not compatible with current implementation.\n";
      abort();
    }
  }
  
  Fopr_Chebyshev* getFoprHerm(const DiracWilsonLike* D){
    kernel_= std::auto_ptr<Fopr_Herm>(ddlinFact_.getFoprHerm(D));
    return new Fopr_Chebyshev(kernel_.get(),N_);
  }
};


#endif
