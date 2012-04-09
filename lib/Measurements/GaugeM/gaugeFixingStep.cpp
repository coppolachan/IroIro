/*! @file gaugeFixStep.cpp */

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

  // even
  GaugeField1D Gtr = max_trace(W_even(Ue,Uo));
  Gtr *= wp;

  for(int site=0; site<Nvh_; ++site){
    SUNmat ut = mat(Gtr,site);
    ut += uwp;
    SetMat(Gtr,ut.reunit(),site);
  }
  gauge_tr_even(Ue,Uo,Gtr);

  // odd
  Gtr = max_trace(W_odd(Ue,Uo));
  Gtr *= wp;

  for(int site=0; site<Nvh_; ++site){
    SUNmat ut = mat(Gtr,site);
    ut += uwp;
    SetMat(Gtr,ut.reunit(),site);
  }
  gauge_tr_odd(Ue,Uo,Gtr);
}

void GaugeFixingStep::gauge_tr_even(GaugeField& Ue,GaugeField& Uo,
				    const GaugeField1D& G)const{
  for(int mu=0; mu<NDIM_; ++mu){

    for(int site=0; site<Nvh_; ++site)
      SetMat(Ue,mat(G,site)*mat(Ue,site,mu),site,mu);

    GaugeField1D Gt = shiftField_oe(G,mu,Backward());

    for(int site=0; site<Nvh_; ++site)
      SetMat(Uo,mat(Uo,site,mu)*mat_dag(Gt,site),site,mu);
  }
}

void GaugeFixingStep::gauge_tr_odd(GaugeField& Ue,GaugeField& Uo,
				   const GaugeField1D& G)const{
  for(int mu=0; mu<NDIM_; ++mu){

    for(int site=0; site<Nvh_; ++site)
      SetMat(Uo,mat(G,site)*mat(Uo,site,mu),site,mu);

    GaugeField1D Gt = shiftField_eo(G,mu,Backward());
    
    for(int site=0; site<Nvh_; ++site)
      SetMat(Ue,mat(Ue,site,mu)*mat_dag(Gt,site),site,mu);
  }
}

const GaugeField1D GaugeFixingStep::W_even(const GaugeField& Ue,
					   const GaugeField& Uo)const{
  GaugeField1D W(Nvh_);
  for(int mu=0; mu<NDIM_-1; ++mu){
    W += DirSlice(Ue,mu);
    GaugeField1D Ut = shiftField_eo(DirSlice(Uo,mu),mu,Forward());

    for(int site=0; site<Nvh_; ++site) AddMat(W,mat_dag(Ut,site),site);
  }
  return W;
}

const GaugeField1D GaugeFixingStep::W_odd(const GaugeField& Ue,
					  const GaugeField& Uo)const{
  GaugeField1D W(Nvh_);
  for(int mu=0; mu<NDIM_-1; ++mu){
    W += DirSlice(Uo,mu);
    GaugeField1D Ut = shiftField_oe(DirSlice(Ue,mu),mu,Forward());

    for(int site=0; site<Nvh_; ++site) AddMat(W,mat_dag(Ut,site),site);
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
}

void GaugeFixingStep::maxTrSU3_1(GaugeField1D& G,GaugeField1D& W)const{

  SUNmat gt;            // will contain SU(2) partial matrix   
  gt.set(2, 0.0,0.0);   //
  gt.set(5, 0.0,0.0);   //      [ w[0] w[1] 0.0 ]
  gt.set(6, 0.0,0.0);   // gt = [ w[3] w[4] 0.0 ]
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
  gt.set(2, 0.0, 0.0);   // gt= [ 0.0  w[4] w[5]]
  gt.set(3, 0.0, 0.0);   //     [ 0.0  w[7] w[8]]
  gt.set(6, 0.0, 0.0);

  for(int site=0; site<W.Nvol(); ++site){

    SUNmat wt = mat(W,site);
    double dst = (wt.r(4) +wt.r(8))*(wt.r(4) +wt.r(8))
                +(wt.i(4) -wt.i(8))*(wt.i(4) -wt.i(8))
                 +(wt.r(7) -wt.r(5))*(wt.r(7) -wt.r(5))
                 +(wt.i(7) +wt.i(5))*(wt.i(7) +wt.i(5));
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
  gt.set(3, 0.0,0.0);   //      [ w[0] 0.0  w[2]]
  gt.set(4, 1.0,0.0);   // gt = [ 0.0  1.0  0.0 ]
  gt.set(5, 0.0,0.0);   //      [ w[6] 0.0  w[8]]
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

