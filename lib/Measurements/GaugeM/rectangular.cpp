//----------------------------------------------------------------------
/*!
  @file rectangular.cpp

  @brief Defines RectangularStaple class functions

*/ 
//----------------------------------------------------------------------

#include "rectangular.hpp"

typedef ShiftField_up<GaugeFieldFormat> FieldUP;
typedef ShiftField_dn<GaugeFieldFormat> FieldDN;
typedef valarray<double> field1d;


Field Staples::upper_1(const Field& g,
		       const Field& cup, 
		       int mu, 
		       int nu) const{
  //       mu,v                               
  //      +-->--+-->--+                                                    
  // nu,w |           |cup_dag(site+mu,nu)
  //  site+     +--<--+                               

  // cup is just a staple, so it is better to call it externally 
  // without the risk or repeating calculations

  valarray<double> w = g[gf_.dir_slice(nu)];
  valarray<double> v = g[gf_.dir_slice(mu)]; 
  FieldUP un(v,sf_,nu);  
  FieldUP cup_mu(cup,sf_,mu);  
  
  valarray<double> res(sf_->size());
  for(int site = 0; site<Nvol_; ++site){
    res[sf_->cslice(0,site)] = (u(w,site)*u(un,site)*u_dag(cup_mu,site)).getva();
  }
  return Field(res);
  
}

Field Staples::upper_1(const ShiftField& gs, 
		       const ShiftField& cup_shift,
		       int mu,
		       int nu) const{
  Field g(gs.getva());
  Field cup(cup_shift.getva());
  return upper_rectangular(g, cup,mu,nu);
}

