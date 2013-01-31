/*!@file fopr_ChebyshevDdagDLin_Factory.hpp
 * @brief factory of fopr_Chebyshev whose kernel is DdagDLin
 This is especially used for the chebyshev acceleration in the eigenmodes calc.
 */
#ifndef FOPRHERMFACTORY_CHEBYSHEVDDAGDLIN_INCLUDED
#define FOPRHERMFACTORY_CHEBYSHEVDDAGDLIN_INCLUDED

#include "include/fopr_Chebyshev.h"
#include <iostream>

class FoprHermFactory_ChebyshevDdagDLin: public FoprHermFactory{
  int N_;
  double slope_;
  double itcpt_;
  std::auto_ptr<Fopr_Herm> ddlin_;  
  /*!<@brief ddlin_ is created, used and destructed within this class.*/
public:
  FoprHermFactory_ChebyshevDdagDLin(const XML::node& node):ddlin_(NULL){

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

      slope_= sgn*2.0/(exu*exu-exl*exl);
      itcpt_= -sgn*(exu*exu+exl*exl)/(exu*exu-exl*exl);
    }else{
      CCIO::cout<<kn_name<<" is not compatible with current implementation.\n";
      abort();
    }
  }
  
  ~FoprHermFactory_ChebyshevDdagDLin(){}

  Fopr_Chebyshev* getFoprHerm(const DiracWilsonLike* D){
    ddlin_= std::auto_ptr<Fopr_Herm>(new Fopr_DdagDLinear(D,slope_,itcpt_)); 

    CCIO::cout<<"Fopr_Chebyshev is going to be constructed\n";
    std::cout.setf(std::ios::showpos);
    CCIO::cout<<"kernel= " <<slope_<<"*DdagD "<<itcpt_<<"\n";
    std::cout.unsetf(std::ios::showpos);

    return new Fopr_Chebyshev(ddlin_.get(),N_);
  }
};

#endif
