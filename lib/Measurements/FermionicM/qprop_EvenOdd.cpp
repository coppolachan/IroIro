//------------------------------------------------------------------------
/*!
 * @file qprop_EvenOdd.cpp
 * @brief Definition of Qprop_EvenOdd member functions 
 */
//-------------------------------------------------------------------------
#include "qprop_EvenOdd.hpp"
#include "Measurements/FermionicM/source.hpp"
#include "Communicator/comm_io.hpp"

void Qprop_EvenOdd::calc(prop_t& xq,Source& src) const{
  xq.clear();
  for(int s=0; s<CommonPrms::instance()->Nd();++s){
    for(int c=0; c<CommonPrms::instance()->Nc();++c){
#if VERBOSITY > 0
      CCIO::cout<<"s="<<s<<" c="<<c<< std::endl;
#endif
      Field sol(fsize_);        
      InvD_.invert(sol,src.mksrc(s,c));
      xq.push_back(sol);
    }
  }
}

void Qprop_EvenOdd::calc(prop_t& xq,Source& src, int s, int c) const{
  Field sol(fsize_);        
  InvD_.invert(sol,src.mksrc(s,c));
  xq.push_back(sol);
}

void Qprop_EvenOdd::calc(prop_t& xq,const prop_t& prp)const{
  xq.clear();

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
#if VERBOSITY > 0
      CCIO::cout<<"s="<<s<<" c="<<c<< std::endl;
#endif
      Field sol(fsize_);        
      InvD_.invert(sol,prp[Nc_*s+c]);
      xq.push_back(sol);
    }
  }
}
