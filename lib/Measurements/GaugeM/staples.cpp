/*! 
  @file staples.cpp
  @brief Defines the staples measurement classes
*/
#include "staples.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include <algorithm>

using namespace FieldUtils;
using namespace SUNmatUtils;

//------------------------------------------------------------
double Staples::plaquette(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaquette called\n");
  return (plaq_s(F) + plaq_t(F))*0.5;
}

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
  return plaq/(Lvol_*NC_*(NDIM_-1));
}

double Staples::plaq_t(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_t called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);

  for(int nu=0; nu < NDIM_-1; ++nu){
    stpl = lower(F,TDIR,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,TDIR)*mat_dag(stpl,site));  // P_zx
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*NC_*(NDIM_-1));
}

double Staples::plaquette_adj(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaquette called\n");
  return (plaq_s_adj(F) + plaq_t_adj(F))*0.5;
}

double Staples::plaq_s_adj(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_s called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);
  SUNmat pl;

  for(int i=0;i<NDIM_-1;++i){
    int j = (i+1)%(NDIM_-1);
    stpl = lower(F,i,j);
    for(int site=0; site<Nvol_; ++site){
      pl = mat(F,site,i)*mat_dag(stpl,site);// P_ij
      double retrace = ReTr(pl);  
      double imtrace = ImTr(pl);  
      plaq += retrace*retrace + imtrace*imtrace - 1.0;
    }
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*(NC_*NC_-1.0)*(NDIM_-1));
}

double Staples::plaq_t_adj(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_t called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);
  SUNmat pl;

  for(int nu=0; nu < NDIM_-1; ++nu){
    stpl = lower(F,TDIR,nu);
    for(int site=0; site<Nvol_; ++site){
      pl = mat(F,site,TDIR)*mat_dag(stpl,site);// P_zx
      double retrace = ReTr(pl);  
      double imtrace = ImTr(pl);  
      plaq += retrace*retrace + imtrace*imtrace - 1.0;
    }
  }
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*(NC_*NC_-1.0)*(NDIM_-1));
}

double Staples::plaq_min(const GaugeField& F,double threshold)const {
  GaugeField1D stpl(Nvol_);

  double pmin(NC_);
  double thldxNc = threshold*NC_;
  int count = 0;
  //double pav =0.0;

  for(int m=0;m<NDIM_;++m){
    for(int n=m+1;n<NDIM_;++n){
      stpl = lower(F,m,n);
      for(int site=0; site<Nvol_; ++site){
	double pl = ReTr(mat(F,site,m)*mat_dag(stpl,site));
	if(pl < thldxNc) count++;

        pmin = std::min(pl,pmin);
	//pav += pl;
      }
    }
  }
  //CCIO::cout<<"averaged plaq= "<<pav/(Lvol_*NC_*6.0)<<"\n";
  int dum;
  int one = Communicator::instance()->reduce_min(pmin,dum,1);

  if(threshold < 1.0){
    double rate(count);
    rate = Communicator::instance()->reduce_sum(rate)/(Lvol_*6.0);
    CCIO::cout<<"   Rate(pl <"<< threshold<<")= "<< rate <<"\n";
  }
  return pmin/NC_;
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
  GaugeField1D c(G.Nvol()); 
  GaugeField1D WupMu = shiftField(w,mu,Forward());
  GaugeField1D VupNu = shiftField(v,nu,Forward());

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
    std::slice isl = c.format.islice(site);
    c.data[isl] = (mat(w,site)*mat(VupNu,site)*mat_dag(WupMu,site)).getva();
    VupNu.data[isl] = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
  }
  c += shiftField(VupNu,nu,Backward());
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

  for(int site=0; site<Nvol_; ++site) 
    c.data[c.format.islice(site)] 
      = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
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

  for(int site=0; site<Nvol_; ++site)
    c.data[c.format.islice(site)] 
      = (mat(w,site)*mat(VupNu,site)*mat_dag(WupMu,site)).getva();
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
//------------------------------------------------------------

GaugeField1D 
Staples::fieldStrength(const GaugeField& G,int mu,int nu) const{
  using namespace Mapping;

  GaugeField1D Vup = upper(G,mu,nu); // Upper staple V_+mu
  GaugeField1D Vdn = lower(G,mu,nu); // Lower staple V_-mu

  GaugeField1D Fmn(Nvol_);
  GaugeField1D Ut(Nvol_);            // temporal objects

  for(int site=0; site<Nvol_; ++site){
    SUNmat u = mat(G,site,mu);    
    SUNmat v = mat_dag(Vup,site); 
    v -= mat_dag(Vdn,site);          // Vup_mu^dag -Vdn_mu^dag
    SetMat(Fmn,(u*v).anti_hermite(),site);
    SetMat(Ut, (v*u).anti_hermite(),site);
  }
  //Fmn +--<--+  Ut +--<--+ 
  //    |     |     |     | 
  // (x)+-->--+     +-->--+(x)
  //    |     |     |     | 
  //    +--<--+     +--<--+  

  Fmn += shiftField(Ut,mu,Backward());
  Fmn *= 0.25;
  return Fmn;
}
