/*! @file gaugeFixStep.cpp 
    @brief defines the gauge fixing step for even/odd sites
*/
#include "gaugeFixingStep.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/macros.hpp"
#include "lib/Main/Geometry/shiftField.hpp"
#include "lib/Main/Geometry/mapping.hpp"
#include  <complex.h>
#ifdef IBM_BGQ_WILSON
#include "Tools/sunAlgUtils.hpp" // this is not really a good option
// please separate the improved routines in a different file
#include <omp.h>
#include "bgqthread.h"
#endif

using namespace SUNmatUtils;
using namespace FieldUtils;
using namespace Mapping;

void GaugeFixingStep::step_SU2(GaugeField& Ue,GaugeField& Uo,
			       OverRelaxing OvR)const{
  GaugeField1D Gt(Nvh_);

  W_even(Gt,Ue,Uo);     /// on the even sites      
  max_trace_ovr(Gt,OvR);
  gauge_tr_even(Ue,Uo,Gt);

  W_odd(Gt,Ue,Uo);      /// on the odd sites      
  max_trace_ovr(Gt,OvR);
  gauge_tr_odd(Ue,Uo,Gt);
}

void GaugeFixingStep::max_trace_ovr(GaugeField1D& G,OverRelaxing OvR)const{

#ifdef IBM_BGQ_WILSON
  double* G_ptr = G.data.getaddr(0); 
#pragma omp parallel 
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC2 = 2*NC_*NC_;

    G_ptr += CC2*ns*omp_get_thread_num();

    for(int i=0; i<ns; ++i){
      max_trace(G_ptr);
      if(OvR) overrelax(G_ptr);
      G_ptr += CC2;
    }
  }
#else
  SUNmat gt;
  for(int site=0; site<Nvh_; ++site){
    max_trace(gt,mat(G,site));
    if(OvR) overrelax(gt);
    SetMat(G,gt,site);
  }
#endif
}

void GaugeFixingStep::overrelax(SUNmat& g)const{
  g *= orp_;
  g += unity()*(1.0-orp_);
  g.reunit();

  /* ///second order overrelaxation
  SUNmat gt = unity()*(1.0-orp_)*(1-orp_);
  gt += g*(2.0*orp_-1.0)*orp_;
  gt += g*g*(orp_-1.0)*orp_;
  g = gt;
  */
}

void GaugeFixingStep::gauge_tr_even(GaugeField& Ue,GaugeField& Uo,
				    const GaugeField1D& Ge)const{
  GaugeField1D  Gp(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* Gp_ptr = Gp.data.getaddr(0);
  double* Ge_ptr = const_cast<GaugeField1D&>(Ge).data.getaddr(0);
  double* Ue_ptr = Ue.data.getaddr(0);
  double* Uo_ptr = Uo.data.getaddr(0);

  GaugeField1D U1d(Nvh_);
  double* U1d_ptr = U1d.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str = CC*ns*omp_get_thread_num();
    const int str2 = str*2;

    for(int mu=0; mu<NDIM_; ++mu){ /*!<gauge tr. is always done for NDIM-dims*/
      shiftField_oe(Gp,Ge_ptr,mu,Forward());

      BGWilsonSU3_MatMult_NN(U1d_ptr+str2,
			     Ge_ptr+str2,Ue_ptr+mu*Nvh_*CC*2+str2,ns);
      BGWilsonSU3_MatEquate(Ue_ptr+mu*Nvh_*CC*2+str2,
			    U1d_ptr+str2,ns);

      BGWilsonSU3_MatMult_ND(U1d_ptr+str2,
			     Uo_ptr+mu*Nvh_*CC*2+str2,Gp_ptr+str2,ns);
      BGWilsonSU3_MatEquate(Uo_ptr+mu*Nvh_*CC*2+str2,
			    U1d_ptr+str2,ns);
    }
  }
#else
  for(int mu=0; mu<NDIM_; ++mu){ /*!<gauge tr. is always done for NDIM-dims*/
    Gp = shiftField_oe(Ge,mu,Forward());
    for(int site=0; site<Nvh_; ++site){
      SetMat(Ue,mat(Ge,site)*mat(Ue,site,mu),site,mu);
      SetMat(Uo,mat(Uo,site,mu)*mat_dag(Gp,site),site,mu);
    }
  }
#endif
}  

void GaugeFixingStep::gauge_tr_odd(GaugeField& Ue,GaugeField& Uo,
				    const GaugeField1D& Go)const{
  GaugeField1D Gp(Nvh_);
  
#ifdef IBM_BGQ_WILSON
  double* Gp_ptr = Gp.data.getaddr(0);
  double* Go_ptr = const_cast<GaugeField1D&>(Go).data.getaddr(0);
  double* Ue_ptr = Ue.data.getaddr(0);
  double* Uo_ptr = Uo.data.getaddr(0);

  GaugeField1D U1d(Nvh_);
  double* U1d_ptr = U1d.data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel 
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str = CC*ns*omp_get_thread_num();
    const int str2 = str*2;

    for(int mu=0; mu<NDIM_; ++mu){ /*!<gauge tr. is always done for NDIM-dims*/
      shiftField_eo(Gp,Go_ptr,mu,Forward());

      BGWilsonSU3_MatMult_ND(U1d_ptr+str2,
			     Ue_ptr+mu*Nvh_*CC*2+str2,Gp_ptr+str2,ns);
      BGWilsonSU3_MatEquate(Ue_ptr+mu*Nvh_*CC*2+str2,
			    U1d_ptr+str2,ns);

      BGWilsonSU3_MatMult_NN(U1d_ptr+str2,
			     Go_ptr+str2,Uo_ptr+mu*Nvh_*CC*2+str2,ns);
      BGWilsonSU3_MatEquate(Uo_ptr+mu*Nvh_*CC*2+str2,
			    U1d_ptr+str2,ns);
    }
  }
#else
  for(int mu=0; mu<NDIM_; ++mu){ /*!<gauge tr. is always done for NDIM-dims*/
    Gp = shiftField_eo(Go,mu,Forward());
    for(int site=0; site<Nvh_; ++site){
      SetMat(Ue,mat(Ue,site,mu)*mat_dag(Gp,site),site,mu);
      SetMat(Uo,mat(Go,site)*mat(Uo,site,mu),site,mu);
    }
  }  
#endif
}

/// this function is for the test purpose
void GaugeFixingStep::gauge_tr(GaugeField& U,const GaugeField1D& G)const{

  for(int mu=0; mu<Ndim_; ++mu){
    const GaugeField1D Gp = shiftField(G,mu,Forward());
    for(int site=0; site<Nvh_; ++site)
      SetMat(U,mat(G,site)*mat(U,site,mu)*mat_dag(Gp,site),site,mu);
  }  
}
//////////////////////

void GaugeFixingStep::step_sdm(GaugeField& Ue,GaugeField& Uo)const{
  GaugeField1D Gtr = gtr_sdm(upu_even(Ue,Uo),umu_even(Ue,Uo));
  gauge_tr_even(Ue,Uo,Gtr);
  Gtr =              gtr_sdm(upu_odd( Ue,Uo),umu_odd( Ue,Uo));
  gauge_tr_odd( Ue,Uo,Gtr);
}

const GaugeField1D GaugeFixingStep::gtr_sdm(const GaugeField1D& upu,
					    const GaugeField1D& umu)const{
  GaugeField1D del = TracelessAntihermite(umu);
  GaugeField1D Gtr(upu.Nvol());

  for(int site=0; site<upu.Nvol(); ++site){
    SUNmat dx = mat(del,site);
    double alpha = ReTr(dx*mat(umu,site));
    alpha /= ReTr(dx*dx*mat(upu,site));
    dx*= alpha*sdmp_;
    SetMat(Gtr,reunit(exponential(dx,12)),site);
  }
  return Gtr;
}

const GaugeField1D GaugeFixingStep::upu_even(const GaugeField& Ue,
					     const GaugeField& Uo)const{
  GaugeField1D upu(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* UpU_ptr = upu.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  GaugeField1D Ute(Nvh_);
  double* Ute_ptr = Ute.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatAdd(UpU_ptr+str2,
			 Ue_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_eo(Ute,Uo_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd(UpU_ptr+str2,
			 Ute_ptr+str2,ns);
    }      
  }
#else
  for(int mu=0; mu<Ndim_; ++mu){
    upu += DirSlice(Ue,mu);
    upu += shiftField_eo(DirSlice(Uo,mu),mu,Backward());
  }
#endif
  return upu;
}

const GaugeField1D GaugeFixingStep::umu_even(const GaugeField& Ue,
					     const GaugeField& Uo)const{
  GaugeField1D umu(Nvh_);
 
#ifdef IBM_BGQ_WILSON
  double* UmU_ptr = umu.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  GaugeField1D Ute(Nvh_);
  double* Ute_ptr = Ute.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatSub(UmU_ptr+str2,
			 Ue_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_eo(Ute,Uo_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd(UmU_ptr+str2,
			 Ute_ptr+str2,ns);
    }      
  }
#else
 for(int mu=0; mu<Ndim_; ++mu){
    umu -= DirSlice(Ue,mu);
    umu += shiftField_eo(DirSlice(Uo,mu),mu,Backward());
  }
#endif
  return umu;
}

const GaugeField1D GaugeFixingStep::upu_odd(const GaugeField& Ue,
					    const GaugeField& Uo)const{
  GaugeField1D upu(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* UpU_ptr = upu.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  GaugeField1D Uto(Nvh_);
  double* Uto_ptr = Uto.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatAdd(UpU_ptr+str2,
			 Uo_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_oe(Uto,Ue_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd(UpU_ptr+str2,
			 Uto_ptr+str2,ns);
    }      
  }
#else
  for(int mu=0; mu<Ndim_; ++mu){
    upu += DirSlice(Uo,mu);
    upu += shiftField_oe(DirSlice(Ue,mu),mu,Backward());
  }
#endif
  return upu;
}

const GaugeField1D GaugeFixingStep::umu_odd(const GaugeField& Ue,
					    const GaugeField& Uo)const{
  GaugeField1D umu(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* UmU_ptr = umu.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  GaugeField1D Uto(Nvh_);
  double* Uto_ptr = Uto.data.getaddr(0);

  BGQThread_Init();
  
#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatSub(UmU_ptr+str2,
			 Uo_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_oe(Uto,Ue_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd(UmU_ptr+str2,
			 Uto_ptr+str2,ns);
    }      
  }
#else
  for(int mu=0; mu<Ndim_; ++mu){
    umu -= DirSlice(Uo,mu);
    umu += shiftField_oe(DirSlice(Ue,mu),mu,Backward());
  }
#endif
  return umu;
}

void GaugeFixingStep::W_even(GaugeField1D& W,
			     const GaugeField& Ue,const GaugeField& Uo)const{
  GaugeField1D Wt(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* W_ptr = W.data.getaddr(0);
  double* Wt_ptr = Wt.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();
    
    BGWilsonSU3_MatZero(W_ptr+str2,ns);

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatAdd(W_ptr+str2,
			 Ue_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_eo(Wt,Uo_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd_ND(W_ptr+str2,
			    Wt_ptr+str2,ns);
    }      
  }
#else
  W = 0.0;
  for(int mu=0; mu<Ndim_; ++mu){
    W += DirSlice(Ue,mu);
    Wt = shiftField_eo(DirSlice(Uo,mu),mu,Backward());
    for(int site=0; site<Nvh_; ++site) AddMat(W,mat_dag(Wt,site),site);
  }
#endif
}

void GaugeFixingStep::W_odd(GaugeField1D& W,
			    const GaugeField& Ue,const GaugeField& Uo)const{
  GaugeField1D Wt(Nvh_);

#ifdef IBM_BGQ_WILSON
  double* W_ptr = W.data.getaddr(0);
  double* Wt_ptr = Wt.data.getaddr(0);
  double* Ue_ptr = const_cast<GaugeField&>(Ue).data.getaddr(0);
  double* Uo_ptr = const_cast<GaugeField&>(Uo).data.getaddr(0);

  BGQThread_Init();

#pragma omp parallel
  {
    const int ns = Nvh_/omp_get_num_threads();
    const int CC = NC_*NC_;
    const int str2 = 2*CC*ns*omp_get_thread_num();
    
    BGWilsonSU3_MatZero(W_ptr+str2,ns);

    for(int mu=0; mu<Ndim_; ++mu){
      BGWilsonSU3_MatAdd(W_ptr+str2,
			 Uo_ptr+mu*Nvh_*CC*2+str2,ns);

      shiftField_oe(Wt,Ue_ptr+mu*Nvh_*CC*2,mu,Backward());

      BGWilsonSU3_MatAdd_ND(W_ptr+str2,
			    Wt_ptr+str2,ns);
    }      
  }
#else
  W = 0.0;
  for(int mu=0; mu<Ndim_; ++mu){
    W += DirSlice(Uo,mu);
    Wt = shiftField_oe(DirSlice(Ue,mu),mu,Backward());
    for(int site=0; site<Nvh_; ++site) AddMat(W,mat_dag(Wt,site),site);
  }
#endif
}

void GaugeFixingStep::max_trace(SUNmat& g,const SUNmat& w)const{
  // Cabibbo-Marinari maximization

  maxTrSU3_1(g,w);
  SUNmat gt;
  maxTrSU3_2(gt,g*w);
  g = gt*g;
  maxTrSU3_3(gt,g*w);
  g = gt*g;
}

void GaugeFixingStep::maxTrSU3_1(SUNmat& gt,const SUNmat& wt)const{
  gt=0.0;         //      [ x[0] x[1] 0.0 ]  
  gt.setr(8,1.0); // gt = [ x[3] x[4] 0.0 ]
                  //      [ 0.0  0.0  1.0 ]
  double dst = (wt.r(0) +wt.r(4))*(wt.r(0) +wt.r(4))
              +(wt.i(0) -wt.i(4))*(wt.i(0) -wt.i(4))
              +(wt.r(1) -wt.r(3))*(wt.r(1) -wt.r(3))
              +(wt.i(1) +wt.i(3))*(wt.i(1) +wt.i(3));
  double fn = 1.0/sqrt(dst);

  gt.set(0, fn*( wt.r(0)+wt.r(4)), fn*(-wt.i(0)+wt.i(4)));
  gt.set(1, fn*(-wt.r(1)+wt.r(3)), fn*(-wt.i(1)-wt.i(3)));
  gt.set(3, fn*( wt.r(1)-wt.r(3)), fn*(-wt.i(1)-wt.i(3)));
  gt.set(4, fn*( wt.r(0)+wt.r(4)), fn*( wt.i(0)-wt.i(4)));
}

void GaugeFixingStep::maxTrSU3_2(SUNmat& gt,const SUNmat& wt)const{
  gt=0.0;         //     [ 1.0  0.0  0.0 ]
  gt.setr(0,1.0); // gt= [ 0.0  x[4] x[5]]       
                  //     [ 0.0  x[7] x[8]]
  double dst = (wt.r(4) +wt.r(8))*(wt.r(4) +wt.r(8))
              +(wt.i(4) -wt.i(8))*(wt.i(4) -wt.i(8))
              +(wt.r(5) -wt.r(7))*(wt.r(5) -wt.r(7))
              +(wt.i(5) +wt.i(7))*(wt.i(5) +wt.i(7));
  double fn = 1.0/sqrt(dst);

  gt.set(4, fn*( wt.r(4)+wt.r(8)), fn*(-wt.i(4)+wt.i(8)));
  gt.set(5, fn*(-wt.r(5)+wt.r(7)), fn*(-wt.i(5)-wt.i(7)));
  gt.set(7, fn*( wt.r(5)-wt.r(7)), fn*(-wt.i(5)-wt.i(7)));
  gt.set(8, fn*( wt.r(4)+wt.r(8)), fn*( wt.i(4)-wt.i(8)));
}

void GaugeFixingStep::maxTrSU3_3(SUNmat& gt,const SUNmat& wt)const{
  gt=0.0;         //      [ x[0] 0.0  x[2]]
  gt.setr(4,1.0); // gt = [ 0.0  1.0  0.0 ]
                  //      [ x[6] 0.0  x[8]]
  double dst = (wt.r(8) +wt.r(0))*(wt.r(8) +wt.r(0))
              +(wt.i(8) -wt.i(0))*(wt.i(8) -wt.i(0))
              +(wt.r(2) -wt.r(6))*(wt.r(2) -wt.r(6))
              +(wt.i(2) +wt.i(6))*(wt.i(2) +wt.i(6));
  double fn = 1.0/sqrt(dst);
    
  gt.set(0, fn*( wt.r(8)+wt.r(0)), fn*( wt.i(8)-wt.i(0)));
  gt.set(2, fn*( wt.r(6)-wt.r(2)), fn*(-wt.i(6)-wt.i(2)));
  gt.set(6, fn*(-wt.r(6)+wt.r(2)), fn*(-wt.i(6)-wt.i(2)));
  gt.set(8, fn*( wt.r(8)+wt.r(0)), fn*(-wt.i(8)+wt.i(0)));
}

/* /// experimental function
void GaugeFixingStep::step_CG(GaugeField& Ue,GaugeField& Uo,
			      GaugeField1D& De,GaugeField1D& Do,
			      GaugeField1D& Pe,GaugeField1D& Po)const{

  GaugeField1D Gtr = gtr_CG(upu_even(Ue,Uo),umu_even(Ue,Uo),De,Pe);
  gauge_tr_even(Ue,Uo,Gtr);
  Gtr =              gtr_CG(upu_odd( Ue,Uo),umu_odd( Ue,Uo),Do,Po);
  gauge_tr_odd( Ue,Uo,Gtr);
}

const GaugeField1D GaugeFixingStep::gtr_CG(const GaugeField1D& upu,
					   const GaugeField1D& umu,
					   GaugeField1D& D,
					   GaugeField1D& P)const{

  GaugeField1D del = TracelessAntihermite(umu);
  GaugeField1D Gtr(upu.Nvol());

  double beta;
  for(int site=0; site<upu.Nvol(); ++site){
    SUNmat dx = mat(del,site);
    beta = ReTr(dx*dx);
    beta /= ReTr(mat(P,site)*mat(D,site).anti_hermite());

    SUNmat px = mat(P,site);
    px *= beta;
    px += dx;
    SetMat(P,px,site);

    double alpha = ReTr(px*mat(umu,site));
    alpha /= ReTr(px*px*mat(upu,site));

    //px*= alpha*sdmp_;
    px*= alpha;      
    SetMat(Gtr,reunit(exponential(px,12)),site);
  }
  D = umu;
  return Gtr;
}
*/

#ifdef IBM_BGQ_WILSON
void GaugeFixingStep::max_trace(double* w)const{

  double _Complex* pW = (double _Complex*)w;
  double _Complex  gt[4],G[9],w1,w2;
  double r1, r2,fn;  

  /*------ maxTrSU3_1 ------*/
  r1 = cabs(conj(*(pW  ))+ *(pW+4));
  r2 = cabs(conj(*(pW+3))- *(pW+1));
  fn =1.0/sqrt(r1*r1 +r2*r2);

  gt[0] = (conj(*(pW  )) + *(pW+4))*fn; //  [ gt[0] gt[1]   0]
  gt[1] = (conj(*(pW+3)) - *(pW+1))*fn; //  [ gt[2] gt[3]   0]
  gt[2] = -conj(gt[1]);                 //  [  0     0    1.0]
  gt[3] =  conj(gt[0]);                 
  /*-----------------------*/

  G[0] = gt[0];  G[1] = gt[1];          /// G = gt*I
  G[3] = gt[2];  G[4] = gt[3];  

  w1 = *(pW  );                         /// w = gt*w
  w2 = *(pW+3);
  *(pW  ) = gt[0]*w1 +gt[1]*w2;
  *(pW+3) = gt[2]*w1 +gt[3]*w2;

  w1 = *(pW+1);
  w2 = *(pW+4);
  *(pW+1) = gt[0]*w1 +gt[1]*w2;
  *(pW+4) = gt[2]*w1 +gt[3]*w2;

  w1 = *(pW+2);
  w2 = *(pW+5);
  *(pW+2) = gt[0]*w1 +gt[1]*w2;
  *(pW+5) = gt[2]*w1 +gt[3]*w2;

  /*------ maxTrSU3_2 ------*/
  r1 = cabs(conj(*(pW+4))+ *(pW+8));
  r2 = cabs(conj(*(pW+7))- *(pW+5));
  fn =1.0/sqrt(r1*r1 +r2*r2);

  gt[0] = (conj(*(pW+4)) + *(pW+8))*fn; //  [1.0    0    0  ]
  gt[1] = (conj(*(pW+7)) - *(pW+5))*fn; //  [  0 gt[0] gt[1]]
  gt[2] = -conj(gt[1]);                 //  [  0 gt[2] gt[3]]
  gt[3] =  conj(gt[0]);                 
  /*-----------------------*/

  G[6] = G[3]*gt[2];                     /// G = gt*G
  G[7] = G[4]*gt[2];
  G[8] = gt[3];
  G[3] *= gt[0];                        
  G[4] *= gt[0];
  G[5] = gt[1];

  w1 = *(pW+3);                          /// w = gt*w
  w2 = *(pW+6);
  *(pW+3) = gt[0]*w1 +gt[1]*w2;
  *(pW+6) = gt[2]*w1 +gt[3]*w2;

  w1 = *(pW+4);
  w2 = *(pW+7);
  *(pW+4) = gt[0]*w1 +gt[1]*w2;
  *(pW+7) = gt[2]*w1 +gt[3]*w2;

  w1 = *(pW+5);
  w2 = *(pW+8);
  *(pW+5) = gt[0]*w1 +gt[1]*w2;
  *(pW+8) = gt[2]*w1 +gt[3]*w2;

  /*------ maxTrSU3_3 ------*/
  r1 = cabs(conj(*(pW  ))+ *(pW+8));
  r2 = cabs(conj(*(pW+6))- *(pW+2));
  fn =1.0/sqrt(r1*r1 +r2*r2);

  gt[0] = (conj(*(pW  )) + *(pW+8))*fn; //  [gt[0]  0  gt[1]]
  gt[1] = (conj(*(pW+6)) - *(pW+2))*fn; //  [  0  1.0   0   ]
  gt[2] = -conj(gt[1]);                 //  [gt[2]  0  gt[3]]  
  gt[3] =  conj(gt[0]);                 
  /*-----------------------*/

  *(pW  ) = gt[0]*G[0] +gt[1]*G[6];    /// w = gt*G : output 
  *(pW+1) = gt[0]*G[1] +gt[1]*G[7];
  *(pW+2) = gt[1]*G[8];
  *(pW+3) = G[3];
  *(pW+4) = G[4];
  *(pW+5) = G[5];
  *(pW+6) = gt[2]*G[0] +gt[3]*G[6];
  *(pW+7) = gt[2]*G[1] +gt[3]*G[7];
  *(pW+8) = gt[3]*G[8];  
}

void GaugeFixingStep::overrelax(double* w)const{

  double _Complex* pW = (double _Complex*)w;
  double _Complex  gt[4],G[9],w1,w2;
  double r1, r2,fn;  

  for(int cc=0; cc<NC_*NC_; ++cc) *(pW+cc) *= orp_;
  for(int c=0; c<NC_; ++c) *(pW+c*(NC_+1)) += 1.0-orp_;

  
  SUNmemAlign::reunit(w);
}
#endif
