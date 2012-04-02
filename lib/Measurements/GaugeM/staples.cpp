/*! 
  @file staples.cpp
  @brief Defines the staples measurement classes
*/
#include "staples.hpp"
#include "Tools/sunMatUtils.hpp"

using namespace FieldUtils;
using namespace SUNmatUtils;

//------------------------------------------------------------
double Staples::plaquette(const GaugeField& F)const {
  CCIO::cout<<"Staples::plaquette called"<<std::endl;
  return (plaq_s(F) + plaq_t(F))*0.5;
}
//------------------------------------------------------------
double Staples::plaq_s(const GaugeField& F)const {
  double plaq = 0.0;
  GaugeField1D stpl;

  for(int i=0;i<NDIM_-1;++i){
    int j = (i+1)%(NDIM_-1);
    stpl = lower(F,i,j);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,i)*mat_dag(stpl,site));  // P_ij
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*NC_*3.0);
}
//------------------------------------------------------------
double Staples::plaq_t(const GaugeField& F)const {
  double plaq = 0.0;
  GaugeField1D stpl;

  for(int nu=0; nu < NDIM_-1; ++nu){
    stpl = lower(F,3,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,TDIR)*mat_dag(stpl,site));  // P_zx
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*NC_*3.0);
}
//------------------------------------------------------------
GaugeField1D Staples::lower(const GaugeField& G, int mu, int nu) const{
  using namespace Mapping;
  //         +     +
  // nu,w_dag|     |w(site+mu,nu) 
  //     site+-->--+ 
  //           mu,v              
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  GaugeField1D c(G.Nvol());
  GaugeField1D WupMu = shiftField(w,mu,Forward());
  for(int site=0; site<Nvol_; ++site) 
    c.data[c.format.islice(site)] 
      = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
    return shiftField(c,nu,Backward());
}
//------------------------------------------------------------
GaugeField1D Staples::upper(const GaugeField& G, int mu, int nu) const{
  using namespace Mapping;
  //       mu,v                               
  //      +-->--+                                                    
  // nu,w |     |t_dag(site+mu,nu)
  //  site+     +                                                          
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  GaugeField1D c(G.Nvol());
  GaugeField1D WupMu = shiftField(w,mu,Forward());
  GaugeField1D VupNu = shiftField(v,nu,Forward());
  for(int site=0; site<Nvol_; ++site)
    c.data[c.format.islice(site)] 
      = (mat(w,site)*mat(VupNu,site)*mat_dag(WupMu,site)).getva();
  return c;
}
//------------------------------------------------------------
void Staples::staple(GaugeField1D& W, const GaugeField& G, int mu) const{
  W = 0.0;
  for(int nu = 0; nu < NDIM_; nu++){
    if(nu != mu){
      W += upper(G,mu,nu);
      W += lower(G,mu,nu);
    }
  }
}

