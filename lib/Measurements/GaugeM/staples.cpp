//----------------------------------------------------------------------
// staples.cpp
//----------------------------------------------------------------------
#include "staples.h"
#include "Tools/sunMatUtils.hpp"

using namespace SUNmat_utils;

typedef ShiftField_up<GaugeFieldFormat> FieldUP;
typedef ShiftField_dn<GaugeFieldFormat> FieldDN;
typedef std::valarray<double> field1d;

//------------------------------------------------------------
double Staples::plaquette(const Field& g)const{ 
  return (plaq_s(g) +plaq_t(g))/2;
}
//------------------------------------------------------------
double Staples::plaquette(const ShiftField& gs)const{
  Field g(gs.getva());
  return (plaq_s(g) +plaq_t(g))/2;
}
//------------------------------------------------------------
double Staples::plaq_s(const Field& g) const{
  double plaq = 0.0;
  //  Field stpl(sf_->size());
  GaugeField1D stpl;

  for(int i=0;i<Ndim_-1;++i){
    int j = (i+1)%(Ndim_-1);
    
    stpl.U = lower(g,i,j);

    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(u(g,gf_,site,i) * u_dag(stpl,site));  // P_ij
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*Nc_*3.0);
}
//------------------------------------------------------------
double Staples::plaq_s(const ShiftField& gs) const{
  Field g(gs.getva());
  return plaq_s(g);
}
//------------------------------------------------------------
double Staples::plaq_t(const Field& g)const{
  double plaq = 0.0;
  GaugeField1D stpl;

  for(int nu=0; nu < Ndim_-1; ++nu){
    stpl.U = lower(g,3,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(u(g,gf_,site,3)*u_dag(stpl,site));  // P_zx
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*Nc_*3.0);
}
//------------------------------------------------------------
double Staples::plaq_t(const ShiftField& gs)const{
  Field g(gs.getva());
  return plaq_t(g);
}
//------------------------------------------------------------
void Staples::staple(Field& W, const Field& g, int mu) const{
  W = 0.0;
  for(int nu = 0; nu < Ndim_; nu++){
    if(nu != mu){
      W += upper(g,mu,nu);
      W += lower(g,mu,nu);
    }
  }
}
//------------------------------------------------------------
void Staples::staple(Field& W, const ShiftField& gs, int mu)const {
  Field g(gs.getva());
  staple(W,g,mu);
}
//------------------------------------------------------------
Field Staples::upper(const Field& g, int mu, int nu) const{

  //       mu,v                               
  //      +-->--+                                                    
  // nu,w |     |t_dag(site+mu,nu)
  //  site+     +                                                             

  field1d w = g[gf_.dir_slice(nu)];
  field1d v = g[gf_.dir_slice(mu)];
  FieldUP um(w,format1d_,mu);
  FieldUP un(v,format1d_,nu);

  field1d c(format1d_->size());
  for(int site=0; site<Nvol_; ++site){
    c[format1d_->cslice(0,site)] = (u(w,*format1d_,site)*u(un,site)*u_dag(um,site)).getva();
  }
  return Field(c);
}
//------------------------------------------------------------		 
Field Staples::upper(const ShiftField& gs, int mu,int nu) const{
  return upper(Field(gs.getva()),mu,nu);
}
//------------------------------------------------------------
Field Staples::lower(const Field& g, int mu, int nu) const{
  //         +     +
  // nu,w_dag|     |w(site+mu,nu) 
  //     site+-->--+ 
  //           mu,v              

  field1d v = g[gf_.dir_slice(mu)];
  field1d w = g[gf_.dir_slice(nu)];

  FieldUP um(w,format1d_,mu);

  field1d c(format1d_->size());
  for(int site=0; site<Nvol_; ++site)
    c[format1d_->cslice(0,site)] = (u_dag(w,*format1d_,site)*u(v,*format1d_,site)*u(um,site)).getva();

  FieldDN un(c,format1d_,nu);

  for(int site=0; site<Nvol_; ++site)
    v[format1d_->cslice(0,site)] = u(un,site).getva();
  return Field(v);
}
//------------------------------------------------------------
Field Staples::lower(const ShiftField& gs, int mu,int nu) const{
  return lower(Field(gs.getva()),mu,nu);
}

