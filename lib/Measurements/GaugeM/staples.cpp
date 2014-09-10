/*! 
  @file staples.cpp
  @brief Defines the staples measurement classes
*/
#include "staples.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"
#include <algorithm>

#ifdef IBM_BGQ_WILSON
#include <omp.h>
#include "bgqthread.h"
#endif

#ifdef SR16K_WILSON
#include "srmwilson.h"
#include "Tools/Architecture_Optimized/srmwilson_cmpl.hpp"
#endif

using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace Mapping;
//------------------------------------------------------------
double Staples::link_trace(const GaugeField& F) const {
  _Message(DEBUG_VERB_LEVEL, "Staples::link_trace called\n");
  double link = 0.0;
  for(int i=0;i<NDIM_;++i){
    for(int site=0; site<Nvol_; ++site)
      link += ReTr(mat(F,site,i));  // P_ij
  }
  link = com_->reduce_sum(link);
  return link/(Lvol_*NC_*NDIM_);
}

double Staples::plaquette(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaquette called\n");
  return (plaq_s(F) + plaq_t(F))*0.5;
}

double Staples::plaq_s(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_s called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);

#ifdef IBM_BGQ_WILSON
  GaugeField1D U_mu, U_nu;
  GaugeField1D UpNu, UpMu;
  GaugeField1D pl;
  //pointers
  double* F_ptr = const_cast<GaugeField&>(F).data.getaddr(0);
  double* stpl_ptr = stpl.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* pl_ptr = pl.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvol_/omp_get_num_threads();
    const int CC2 = NC_*NC_*2;
    const int str2 = CC2*ns*omp_get_thread_num();

    for(int mu=0; mu<NDIM_-1; ++mu){
      BGWilsonSU3_MatEquate(U_mu_ptr+str2,F_ptr+mu*Nvol_*CC2+str2,ns); //U_mu links
      int nu = (mu+1)%(NDIM_-1);
      BGWilsonSU3_MatEquate(U_nu_ptr+str2,F_ptr+nu*Nvol_*CC2+str2,ns); //U_nu links

      shiftField(UpMu,U_nu_ptr,mu,Forward());
      shiftField(UpNu,U_mu_ptr,nu,Forward());
      // upper staple
      BGWilsonSU3_MatMult_NND(stpl_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
      // plaquette 
      BGWilsonSU3_MatMult_ND(pl_ptr+str2,U_mu_ptr+str2,stpl_ptr+str2,ns);
#pragma omp for reduction(+:plaq)
      for(int site=0; site<Nvol_; ++site)
	  plaq += pl_ptr[CC2*site]+pl_ptr[8+CC2*site]+pl_ptr[16+CC2*site]; ///ReTr
    }
  }
#else  
#ifdef SR16K_WILSON
  GaugeField1D U_mu, U_nu;
  GaugeField1D UpNu, UpMu;
  GaugeField1D pl;
  //pointers
  double* F_ptr = const_cast<GaugeField&>(F).data.getaddr(0);
  double* stpl_ptr = stpl.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* pl_ptr = pl.data.getaddr(0);

  const int CC2 = NC_*NC_*2;

  for(int mu=0; mu<NDIM_-1; ++mu){
    SRCMPL::SRWilsonSU3_MatEquate(U_mu_ptr,F_ptr+mu*Nvol_*CC2,Nvol_);//U_mu links
    int nu = (mu+1)%(NDIM_-1);
    SRCMPL::SRWilsonSU3_MatEquate(U_nu_ptr,F_ptr+nu*Nvol_*CC2,Nvol_);//U_nu links

    shiftField(UpMu,U_nu_ptr,mu,Forward());
    shiftField(UpNu,U_mu_ptr,nu,Forward());
    // upper staple
    SRWilsonSU3_MatMult_NND(stpl_ptr,U_nu_ptr,UpNu_ptr,UpMu_ptr,Nvol_);
    // plaquette 
    SRWilsonSU3_MatMult_ND(pl_ptr,U_mu_ptr,stpl_ptr,Nvol_);

    for(int site=0; site<Nvol_; ++site)
      plaq += pl_ptr[CC2*site]+pl_ptr[8+CC2*site]+pl_ptr[16+CC2*site]; ///ReTr
  }
#else
  for(int i=0;i<NDIM_-1;++i){
    int j = (i+1)%(NDIM_-1);
    stpl = lower(F,i,j);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,i)*mat_dag(stpl,site));  // P_ij
  }
#endif
#endif
  plaq = com_->reduce_sum(plaq);
  return plaq/(Lvol_*NC_*(NDIM_-1));
}

double Staples::plaq_t(const GaugeField& F)const {
  _Message(DEBUG_VERB_LEVEL, "Staples::plaq_t called\n");
  double plaq = 0.0;
  GaugeField1D stpl(Nvol_);

#ifdef IBM_BGQ_WILSON
  GaugeField1D U_mu, U_nu;
  GaugeField1D UpNu, UpMu;
  GaugeField1D pl;
  //pointers
  double* F_ptr = const_cast<GaugeField&>(F).data.getaddr(0);
  double* stpl_ptr = stpl.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* pl_ptr = pl.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvol_/omp_get_num_threads();
    const int CC2 = NC_*NC_*2;
    const int str2 = CC2*ns*omp_get_thread_num();

    int mu=TDIR;
    BGWilsonSU3_MatEquate(U_mu_ptr+str2,F_ptr+mu*Nvol_*CC2+str2,ns); //U_mu links
    for(int nu=0; nu<NDIM_-1; ++nu){
      BGWilsonSU3_MatEquate(U_nu_ptr+str2,F_ptr+nu*Nvol_*CC2+str2,ns); //U_nu links

      shiftField(UpMu,U_nu_ptr,mu,Forward());
      shiftField(UpNu,U_mu_ptr,nu,Forward());
      // upper staple
      BGWilsonSU3_MatMult_NND(stpl_ptr+str2,U_nu_ptr+str2,UpNu_ptr+str2,UpMu_ptr+str2,ns);
      // plaquette 
      BGWilsonSU3_MatMult_ND(pl_ptr+str2,U_mu_ptr+str2,stpl_ptr+str2,ns);
#pragma omp for reduction(+:plaq)
      for(int site=0; site<Nvol_; ++site)
	  plaq += pl_ptr[CC2*site]+pl_ptr[8+CC2*site]+pl_ptr[16+CC2*site]; ///ReTr
    }
  }
#else
#ifdef SR16K_WILSON
  GaugeField1D U_mu, U_nu;
  GaugeField1D UpNu, UpMu;
  GaugeField1D pl;
  //pointers
  double* F_ptr = const_cast<GaugeField&>(F).data.getaddr(0);
  double* stpl_ptr = stpl.data.getaddr(0);
  double* U_mu_ptr = U_mu.data.getaddr(0);
  double* U_nu_ptr = U_nu.data.getaddr(0);
  double* UpMu_ptr = UpMu.data.getaddr(0);
  double* UpNu_ptr = UpNu.data.getaddr(0);
  double* pl_ptr = pl.data.getaddr(0);

  const int CC2 = NC_*NC_*2;

  int mu=TDIR;
  SRCMPL::SRWilsonSU3_MatEquate(U_mu_ptr,F_ptr+mu*Nvol_*CC2,Nvol_); //U_mu links
  for(int nu=0; nu<NDIM_-1; ++nu){
    SRCMPL::SRWilsonSU3_MatEquate(U_nu_ptr,F_ptr+nu*Nvol_*CC2,Nvol_); //U_nu links

    shiftField(UpMu,U_nu_ptr,mu,Forward());
    shiftField(UpNu,U_mu_ptr,nu,Forward());
    // upper staple
    SRWilsonSU3_MatMult_NND(stpl_ptr,U_nu_ptr,UpNu_ptr,UpMu_ptr,Nvol_);
    // plaquette 
    SRWilsonSU3_MatMult_ND(pl_ptr,U_mu_ptr,stpl_ptr,Nvol_);

    for(int site=0; site<Nvol_; ++site)
      plaq += pl_ptr[CC2*site]+pl_ptr[8+CC2*site]+pl_ptr[16+CC2*site]; ///ReTr
  }
#else
  for(int nu=0; nu < NDIM_-1; ++nu){
    stpl = lower(F,TDIR,nu);
    for(int site=0; site<Nvol_; ++site)
      plaq += ReTr(mat(F,site,TDIR)*mat_dag(stpl,site));  // P_zx
  }
#endif
#endif
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

GaugeField1D Staples::upper_lower(const GaugeField& G,int mu,int nu,const Field aux)const{
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
  SUNmat temp;

  int mu_aux = mu;
  int nu_aux = nu;
  int limit;
  double sign = 1.0;
  int sector = -1;
  
  if (mu>nu){
    mu_aux = nu;
    nu_aux = mu;
    sign = -1.0;
  }
  for(int m=0; m<= mu_aux; m++){
    if(m<mu_aux) limit = NDIM_-1;
    else         limit = nu_aux;
    for(int n=m+1; n<=limit; n++) sector++;
  }
  for(int site=0; site<Nvol_; ++site){
    std::slice isl = c.format.islice(site);
    //std::complex<double>     fact(aux[2*site+2*Nvol_*sector],  sign*aux[2*site+1+2*Nvol_*sector]);
    //std::complex<double> fact_dag(aux[2*site+2*Nvol_*sector], -sign*aux[2*site+1+2*Nvol_*sector]);
    temp = mat(w,site)* aux[2*site+2*Nvol_*sector];
    c.data[isl] = (temp*mat(VupNu,site)*mat_dag(WupMu,site)).getva();    
    temp = mat_dag(w,site)*aux[2*site+2*Nvol_*sector];
    VupNu.data[isl] = (temp*mat(v,site)*mat(WupMu,site)).getva();
  }
  c += shiftField(VupNu,nu,Backward());

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
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* w_ptr = w.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvol_/omp_get_num_threads();
    const int CC2 = NC_*NC_*2;
    const int str2 = CC2*ns*omp_get_thread_num();

    BGWilsonSU3_MatEquate(v_ptr+str2,G_ptr+mu*Nvol_*CC2+str2,ns); //U_mu links
    BGWilsonSU3_MatEquate(w_ptr+str2,G_ptr+nu*Nvol_*CC2+str2,ns); //U_nu links
    shiftField(WupMu,w_ptr,mu,Forward());
    BGWilsonSU3_MatMult_DNN(c_ptr+str2,w_ptr+str2,v_ptr+str2,WupMu_ptr+str2,ns);
  }
#else
#ifdef SR16K_WILSON 
  GaugeField1D v,w,c, WupMu;
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);
  double* w_ptr = w.data.getaddr(0);

  const int CC2 = NC_*NC_*2;
  
  SRCMPL::SRWilsonSU3_MatEquate(v_ptr,G_ptr+mu*Nvol_*CC2,Nvol_); //U_mu links
  SRCMPL::SRWilsonSU3_MatEquate(w_ptr,G_ptr+nu*Nvol_*CC2,Nvol_); //U_nu links
  shiftField(WupMu,w_ptr,mu,Forward());
  SRWilsonSU3_MatMult_DNN(c_ptr,w_ptr,v_ptr,WupMu_ptr,Nvol_);
#else
  GaugeField1D v = DirSlice(G,mu);
  GaugeField1D w = DirSlice(G,nu);
  GaugeField1D c(G.Nvol());
  GaugeField1D WupMu = shiftField(w,mu,Forward());

  for(int site=0; site<Nvol_; ++site) 
    c.data[c.format.islice(site)] 
      = (mat_dag(w,site)*mat(v,site)*mat(WupMu,site)).getva();
#endif
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
  GaugeField1D c,v,WupMu,VupNu;
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvol_/omp_get_num_threads();
    const int CC2 = NC_*NC_*2;
    const int str2 = CC2*ns*omp_get_thread_num();

    BGWilsonSU3_MatEquate(v_ptr+str2,G_ptr+mu*Nvol_*CC2+str2,ns); //U_mu links
    shiftField(VupNu,v_ptr,nu,Forward());

    BGWilsonSU3_MatEquate(v_ptr+str2,G_ptr+nu*Nvol_*CC2+str2,ns); //U_nu links
    shiftField(WupMu,v_ptr,mu,Forward());

    BGWilsonSU3_MatMult_NND(c_ptr+str2,v_ptr+str2,VupNu_ptr+str2,WupMu_ptr+str2,ns);
  }
#else
#ifdef SR16K_WILSON 
  GaugeField1D c,v,WupMu,VupNu;
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* c_ptr = c.data.getaddr(0);
  double* v_ptr = v.data.getaddr(0);
  double* VupNu_ptr = VupNu.data.getaddr(0);
  double* WupMu_ptr = WupMu.data.getaddr(0);

  const int CC2 = NC_*NC_*2;
  
  SRCMPL::SRWilsonSU3_MatEquate(v_ptr,G_ptr+mu*Nvol_*CC2,Nvol_); //U_mu links
  shiftField(VupNu,v_ptr,nu,Forward());

  SRCMPL::SRWilsonSU3_MatEquate(v_ptr,G_ptr+nu*Nvol_*CC2,Nvol_); //U_nu links
  shiftField(WupMu,v_ptr,mu,Forward());

  SRWilsonSU3_MatMult_NND(c_ptr,v_ptr,VupNu_ptr,WupMu_ptr,Nvol_);
  
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
  
#ifdef IBM_BGQ_WILSON
  // temporal fields
  GaugeField1D U1;
  // pointers
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* Vup_ptr = Vup.data.getaddr(0);
  double* Vdn_ptr = Vdn.data.getaddr(0);
  double* Fmn_ptr = Fmn.data.getaddr(0);
  double* Ut_ptr = Ut.data.getaddr(0);

  double* U1_ptr = U1.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {  
    const int ns = Nvol_/omp_get_num_threads();
    const int CC2 = NC_*NC_*2;
    const int str2 = CC2*ns*omp_get_thread_num();
    
    BGWilsonSU3_MatEquate(U1_ptr+str2,G_ptr+mu*Nvol_*CC2+str2,ns);//U_mu links 
    BGWilsonSU3_MatSub(Vup_ptr+str2,Vdn_ptr+str2,ns);             //Left leaves

    BGWilsonSU3_MatMult_ND(Fmn_ptr+str2,U1_ptr+str2,Vup_ptr+str2,ns);
    BGWilsonSU3_MatMult_DN(Ut_ptr+str2,Vup_ptr+str2,U1_ptr+str2,ns);
    shiftField(U1, Ut_ptr,mu,Backward());

    BGWilsonSU3_MatAdd(Fmn_ptr+str2,U1_ptr+str2,ns); // Fmn 
  }
#else
#ifdef SR16K_WILSON
  // temporal fields
  GaugeField1D U1;
  // pointers
  double* G_ptr = const_cast<GaugeField&>(G).data.getaddr(0);
  double* Vup_ptr = Vup.data.getaddr(0);
  double* Vdn_ptr = Vdn.data.getaddr(0);
  double* Fmn_ptr = Fmn.data.getaddr(0);
  double* Ut_ptr = Ut.data.getaddr(0);
  double* U1_ptr = U1.data.getaddr(0);
  const int CC2 = NC_*NC_*2;

  SRCMPL::SRWilsonSU3_MatEquate(U1_ptr,G_ptr+mu*Nvol_*CC2,Nvol_); // U_mu links 
  SRCMPL::SRWilsonSU3_MatSub(Vup_ptr,Vdn_ptr,Nvol_); // Left leaves

  SRWilsonSU3_MatMult_ND(Fmn_ptr,U1_ptr, Vup_ptr,Nvol_);
  SRWilsonSU3_MatMult_DN(Ut_ptr, Vup_ptr,U1_ptr, Nvol_);
  shiftField(U1,Ut_ptr,mu,Backward());
  
  SRCMPL::SRWilsonSU3_MatAdd(Fmn_ptr,U1_ptr, Nvol_); // Fmn 
#else
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

#endif
#endif
  Fmn = TracelessAntihermite(Fmn);  // traceless, anti-hermite
  Fmn *= 0.25;
  return Fmn;
}
