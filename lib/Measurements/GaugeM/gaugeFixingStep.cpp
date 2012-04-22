/*! @file gaugeFixStep.cpp 
    @brief defines the gauge fixing step for even/odd sites
*/
#include "gaugeFixingStep.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/macros.hpp"
#include "lib/Main/Geometry/mapping.hpp"
using namespace SUNmatUtils;
using namespace FieldUtils;
using namespace Mapping;

void GaugeFixingStep::gfix_step(GaugeField& Ue,GaugeField& Uo,
				double wp)const{
  SUNmat uwp = unity();
  uwp *= 1.0-wp;

 // on the even sites
  GaugeField1D Gt = max_trace(W_even(Ue,Uo));
  Gt *= wp;  
  for(int site=0; site<Gt.Nvol(); ++site)
    SetMat(Gt,reunit(mat(Gt,site) +uwp),site);
  gauge_tr_even(Ue,Uo,Gt);

  // on the odd sites
  Gt = max_trace(W_odd(Ue,Uo));
  Gt *= wp;
  for(int site=0; site<Gt.Nvol(); ++site)
    SetMat(Gt,reunit(mat(Gt,site) +uwp),site);
  gauge_tr_odd(Ue,Uo,Gt);
}

// this function is just for test
void GaugeFixingStep::gauge_tr(GaugeField& U,const GaugeField1D& G)const{

  for(int mu=0; mu<NDIM_; ++mu){
    const GaugeField1D Gp = shiftField(G,mu,Forward());
    for(int site=0; site<G.Nvol(); ++site)
      SetMat(U,mat(G,site)*mat(U,site,mu)*mat_dag(Gp,site),site,mu);
  }
}  

void GaugeFixingStep::gauge_tr_even(GaugeField& Ue,GaugeField& Uo,
				    const GaugeField1D& Ge)const{
  for(int mu=0; mu<NDIM_; ++mu){
    GaugeField1D Gp = shiftField_oe(Ge,mu,Forward());
    for(int site=0; site<Ge.Nvol(); ++site){
      SetMat(Ue,mat(Ge,site)*mat(Ue,site,mu),site,mu);
      SetMat(Uo,mat(Uo,site,mu)*mat_dag(Gp,site),site,mu);
    }
  }
}  

void GaugeFixingStep::gauge_tr_odd(GaugeField& Ue,GaugeField& Uo,
				    const GaugeField1D& Go)const{
  for(int mu=0; mu<NDIM_; ++mu){
    GaugeField1D Gp = shiftField_eo(Go,mu,Forward());
    for(int site=0; site<Go.Nvol(); ++site){
      SetMat(Ue,mat(Ue,site,mu)*mat_dag(Gp,site),site,mu);
      SetMat(Uo,mat(Go,site)*mat(Uo,site,mu),site,mu);
    }
  }  
}

const GaugeField1D GaugeFixingStep::W_even(const GaugeField& Ue,
					   const GaugeField& Uo)const{
  GaugeField1D W(Uo.Nvol());
  for(int mu=0; mu<Ndim_; ++mu){
    W += DirSlice(Ue,mu);
    GaugeField1D Wt = shiftField_eo(DirSlice(Uo,mu),mu,Backward());
    for(int site=0; site<Uo.Nvol(); ++site) AddMat(W,mat_dag(Wt,site),site);
  }
  return W;
}

const GaugeField1D GaugeFixingStep::W_odd(const GaugeField& Ue,
					  const GaugeField& Uo)const{
  GaugeField1D W(Ue.Nvol());
  for(int mu=0; mu<Ndim_; ++mu){
    W += DirSlice(Uo,mu);
    GaugeField1D Wt = shiftField_oe(DirSlice(Ue,mu),mu,Backward());
    for(int site=0; site<Uo.Nvol(); ++site) AddMat(W,mat_dag(Wt,site),site);
  }
  return W;
}

const GaugeField1D GaugeFixingStep::max_trace(const GaugeField1D& W)const{
  // Cabibbo-Marinari maximization
  SUNmat uni = unity();
  GaugeField1D G(W.Nvol());  
  
  for(int site=0; site<W.Nvol(); ++site) SetMat(G,uni,site);

  GaugeField1D Wt = W;
  assert(NC_==3); // assuming SU(3) gauge theory only

  maxTrSU3_1(G,Wt);
  maxTrSU3_2(G,Wt);
  maxTrSU3_3(G,Wt);
  return G;
}

void GaugeFixingStep::maxTrSU3_1(GaugeField1D& G,GaugeField1D& W)const{

  SUNmat gt;            // will contain SU(2) partial matrix   
  gt.set(2, 0.0,0.0);   //
  gt.set(5, 0.0,0.0);   //      [ x[0] x[1] 0.0 ]
  gt.set(6, 0.0,0.0);   // gt = [ x[3] x[4] 0.0 ]
  gt.set(7, 0.0,0.0);   //      [ 0.0  0.0  1.0 ]
  gt.set(8, 1.0,0.0);
    
  for(int site=0; site<W.Nvol(); ++site){

    SUNmat wt = mat(W,site);
    double dst = (wt.r(0) +wt.r(4))*(wt.r(0) +wt.r(4))
                +(wt.i(0) -wt.i(4))*(wt.i(0) -wt.i(4))
                +(wt.r(1) -wt.r(3))*(wt.r(1) -wt.r(3))
                +(wt.i(1) +wt.i(3))*(wt.i(1) +wt.i(3));
    double fn = 1.0/sqrt(dst);

    gt.set(0, fn*( wt.r(0)+wt.r(4)), fn*(-wt.i(0)+wt.i(4)));
    gt.set(1, fn*(-wt.r(1)+wt.r(3)), fn*(-wt.i(1)-wt.i(3)));
    gt.set(3, fn*( wt.r(1)-wt.r(3)), fn*(-wt.i(1)-wt.i(3)));
    gt.set(4, fn*( wt.r(0)+wt.r(4)), fn*( wt.i(0)-wt.i(4)));

    SetMat(W,gt*wt,site);
    SetMat(G,gt*mat(G,site),site);
  }
}

void GaugeFixingStep::maxTrSU3_2(GaugeField1D& G,GaugeField1D& W)const{

  SUNmat gt;             // will contain SU(2) partial matrix   
  gt.set(0, 1.0, 0.0);   //
  gt.set(1, 0.0, 0.0);   //     [ 1.0  0.0  0.0 ]
  gt.set(2, 0.0, 0.0);   // gt= [ 0.0  x[4] x[5]]
  gt.set(3, 0.0, 0.0);   //     [ 0.0  x[7] x[8]]
  gt.set(6, 0.0, 0.0);

  for(int site=0; site<W.Nvol(); ++site){

    SUNmat wt = mat(W,site);
    double dst = (wt.r(4) +wt.r(8))*(wt.r(4) +wt.r(8))
                +(wt.i(4) -wt.i(8))*(wt.i(4) -wt.i(8))
                +(wt.r(5) -wt.r(7))*(wt.r(5) -wt.r(7))
                +(wt.i(5) +wt.i(7))*(wt.i(5) +wt.i(7));
    double fn = 1.0/sqrt(dst);

    gt.set(4, fn*( wt.r(4)+wt.r(8)), fn*(-wt.i(4)+wt.i(8)));
    gt.set(5, fn*(-wt.r(5)+wt.r(7)), fn*(-wt.i(5)-wt.i(7)));
    gt.set(7, fn*( wt.r(5)-wt.r(7)), fn*(-wt.i(5)-wt.i(7)));
    gt.set(8, fn*( wt.r(4)+wt.r(8)), fn*( wt.i(4)-wt.i(8)));

    SetMat(W,gt*wt,site);
    SetMat(G,gt*mat(G,site),site);
  }
}

void GaugeFixingStep::maxTrSU3_3(GaugeField1D& G,GaugeField1D& W)const{

  SUNmat gt;            // will contain SU(2) partial matrix   
  gt.set(1, 0.0,0.0);   //
  gt.set(3, 0.0,0.0);   //      [ x[0] 0.0  x[2]]
  gt.set(4, 1.0,0.0);   // gt = [ 0.0  1.0  0.0 ]
  gt.set(5, 0.0,0.0);   //      [ x[6] 0.0  x[8]]
  gt.set(7, 0.0,0.0); 

  for(int site=0; site<W.Nvol(); ++site){

    SUNmat wt = mat(W,site);
    double dst = (wt.r(8) +wt.r(0))*(wt.r(8) +wt.r(0))
                +(wt.i(8) -wt.i(0))*(wt.i(8) -wt.i(0))
                +(wt.r(2) -wt.r(6))*(wt.r(2) -wt.r(6))
                +(wt.i(2) +wt.i(6))*(wt.i(2) +wt.i(6));
    double fn = 1.0/sqrt(dst);
    
    gt.set(0, fn*( wt.r(8)+wt.r(0)), fn*( wt.i(8)-wt.i(0)));
    gt.set(2, fn*( wt.r(6)-wt.r(2)), fn*(-wt.i(6)-wt.i(2)));
    gt.set(6, fn*(-wt.r(6)+wt.r(2)), fn*(-wt.i(6)-wt.i(2)));
    gt.set(8, fn*( wt.r(8)+wt.r(0)), fn*(-wt.i(8)+wt.i(0)));

    SetMat(W,gt*wt,site);
    SetMat(G,gt*mat(G,site),site);
  }
}

