/*!
 * @file qprop_optimalDomainWall.cpp
 *
 * @brief Declaration of Qprop optimal domain wall fermions QpropDWF
 *
 */

#include "qprop_DomainWall.hpp"
#include "Communicator/comm_io.hpp"

void QpropDWF::calc(prop_t& xq,Source& src) const{
  xq.clear();

  for(int s=0; s<Nd_;++s) {
     for(int c=0; c<Nc_;++c) {
       CCIO::cout << "Dirac = " << s << " Color = " << c << std::endl;
       xq.push_back(Dgw_.mult_inv(src.mksrc(s,c)));
     }
  }
}

//for testing purposes
void QpropDWF::calc(prop_t& xq,
		    Source& src,
		    int dirac_index,
		    int color) const{
  xq.clear();
  CCIO::cout << "Dirac = " << dirac_index 
	    << " Color = " << color << std::endl;
  xq.push_back(Dgw_.mult_inv(src.mksrc(dirac_index,color)));

}

void QpropDWF::calc(prop_t& xq,const prop_t& prp)const{
  xq.clear();
  
  for(int s=0; s<Nd_;++s) {
    for(int c=0; c<Nc_;++c) {
      CCIO::cout << "Dirac = " << s << " Color = " << c << std::endl;
      xq.push_back(Dgw_.mult_inv(prp[Nc_*s+c]));
    }
  }
}
