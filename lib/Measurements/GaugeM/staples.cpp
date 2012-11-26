/*! 
  @file staples.cpp
  @brief Defines the staples measurement classes
*/
#include "staples.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

using namespace FieldUtils;
using namespace SUNmatUtils;



//------------------------------------------------------------
double Staples::plaquette(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaquette called\n");
  return (plaq_s(F) + plaq_t(F))*0.5;
}
//------------------------------------------------------------
double Staples::plaq_s(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_s called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);

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
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_t called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);

  for(int nu=0; nu < NDIM_-1; ++nu){
    stpl = lower(F,3,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,TDIR)*mat_dag(stpl,site));  // P_zx
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*NC_*3.0);
}
//------------------------------------------------------------
GaugeField1D Staples::upper_lower(const GaugeField& G, int mu, int nu) const{
  _Message(DEBUG_VERB_LEVEL, "Staples::upper_lower called\n");
  using namespace Mapping;
  //       mu,v                               
  //      +-->--+                                                    
  // nu,w |     |t_dag(site+mu,nu)
  //  site+     +
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  CCIO::cout << "w*v  = " << w.data*v.data << "\n";
  GaugeField1D c(G.Nvol()); 
  GaugeField1D WupMu = shiftField(w,mu,Forward());
  GaugeField1D VupNu = shiftField(v,nu,Forward());

  CCIO::cout << "WupMu*v  = " << WupMu.data*v.data << "\n";
  CCIO::cout << "VupNu*v  = " << VupNu.data*v.data << "\n";
  CCIO::cout << "VupNu*WupMu  = " << VupNu.data*WupMu.data << "\n";

#ifdef IBM_BGQ_WILSON   
  double* c_ptr = c.data.getaddr(0);
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* w_ptr = w.data.getaddr(0);

  BGWilsonSU3_MatMult_NND(c_ptr, w_ptr, VupNu_ptr, WupMu_ptr, Nvol_);
  BGWilsonSU3_MatMult_DNN(VupNu_ptr, w_ptr, v_ptr, WupMu_ptr, Nvol_);
  c += shiftField(VupNu,nu,Backward());
#else
  for(int site=0; site<Nvol_; ++site){
    c.data[c.format.islice(site)] = (mat(w,site)*mat(VupNu,site)*mat_dag(WupMu,site)).getva();
  }
  CCIO::cout << "c norm: " << c.norm() << "\n";

  for(int site=0; site<Nvol_; ++site) {
    VupNu.data[c.format.islice(site)] = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
  }

  GaugeField1D temp = shiftField(VupNu,nu,Backward());
  CCIO::cout << "temp*VupNu = " << temp.data*VupNu.data << "\n";
  CCIO::cout << "temp*c     = " << temp.data*c.data << "\n";
  CCIO::cout << "VupNu*c    = " << VupNu.data*c.data << "\n";
  c += temp;
  CCIO::cout << "2- c norm: " << c.norm() << "mu,nu  " << mu << " " << nu << "\n";
#endif
  return c;
}


GaugeField1D Staples::lower(const GaugeField& G, int mu, int nu) const{
  _Message(DEBUG_VERB_LEVEL, "Staples::lower called\n");
  using namespace Mapping;
  //         +     +
  // nu,w_dag|     |w(site+mu,nu) 
  //     site+-->--+ 
  //           mu,v              

#ifdef IBM_BGQ_WILSON 
  GaugeField1D v,w,c, WupMu;
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* w_ptr = w.data.getaddr(0);

  DirSliceBGQ(v,G,mu);
  DirSliceBGQ(w,G,nu);
  shiftField(WupMu,w_ptr,mu,Forward());


  BGWilsonSU3_MatMult_DNN(c_ptr, w_ptr, v_ptr, WupMu_ptr, Nvol_);
#else
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  GaugeField1D c(G.Nvol());
  GaugeField1D WupMu = shiftField(w,mu,Forward());



  for(int site=0; site<Nvol_; ++site) {
    c.data[c.format.islice(site)] = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
  }


   #endif
 

  return shiftField(c,nu,Backward());
}
//------------------------------------------------------------
GaugeField1D Staples::upper(const GaugeField& G, int mu, int nu) const{
  _Message(DEBUG_VERB_LEVEL, "Staples::upper called\n");
  using namespace Mapping;
  //       mu,v                               
  //      +-->--+                                                    
  // nu,w |     |t_dag(site+mu,nu)
  //  site+     +              
 
#ifdef IBM_BGQ_WILSON 
  GaugeField1D v, WupMu, VupNu;
  GaugeField1D c(G.Nvol());
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);

  DirSliceBGQ(v,G,mu); 
  shiftField(VupNu,v_ptr,nu,Forward());
  DirSliceBGQ(v,G,nu);
  shiftField(WupMu,v_ptr,mu,Forward());
 
  BGWilsonSU3_MatMult_NND(c_ptr, v_ptr, VupNu_ptr, WupMu_ptr, Nvol_);
#else
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  GaugeField1D c(G.Nvol());
  GaugeField1D WupMu = shiftField(w,mu,Forward());
  GaugeField1D VupNu = shiftField(v,nu,Forward());

  for(int site=0; site<Nvol_; ++site){
    c.data[c.format.islice(site)] = (mat(w,site)*mat(VupNu,site)*mat_dag(WupMu,site)).getva();
  }

   #endif

  return c;
}
//------------------------------------------------------------
void Staples::staple(GaugeField1D& W, const GaugeField& G, int mu) const{
  _Message(DEBUG_VERB_LEVEL, "Staples::staple called\n");
  W = 0.0;
  for(int nu = 0; nu < NDIM_; nu++){
    if(nu != mu){
      W += upper(G,mu,nu);
      W += lower(G,mu,nu);
    }
  }
}


