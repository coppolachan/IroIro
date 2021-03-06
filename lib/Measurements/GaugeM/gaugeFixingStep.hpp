/*! @file gaugeFixingStep.hpp
    @brief this class defines one step of the gauge fixing, which is 
    common to Coulomb and Landau gauge.
 */
#ifndef GAUGEFIXINGSTEP_INCLUDED
#define GAUGEFIXINGSTEP_INCLUDED

#include "include/common_fields.hpp"
#include "include/macros.hpp"
#include "Tools/sunMat.hpp"

enum Condition{Coulomb,Landau};
enum OverRelaxing{ORoff,ORon};

class GaugeFixingStep{
private:
  int Nvh_;  
  int Ndim_;  
  double orp_,sdmp_;

  void W_even(GaugeField1D& Gt,const GaugeField& Ue,const GaugeField& Uo)const;
  void W_odd( GaugeField1D& Gt,const GaugeField& Ue,const GaugeField& Uo)const;

  const GaugeField1D upu_even(const GaugeField& Ue,const GaugeField& Uo)const;
  const GaugeField1D upu_odd( const GaugeField& Ue,const GaugeField& Uo)const;
  const GaugeField1D gtr_sdm(const GaugeField1D& upu,const GaugeField1D& umu)const;  

#if defined IBM_BGQ_WILSON || defined SR16K_WILSON
  void max_trace(double* w)const;
  void overrelax(double* g)const;
#else
  void max_trace(SUNmat& g,const SUNmat& w)const;
  void overrelax(SUNmat& g)const;

  void maxTrSU3_1(SUNmat& gt,const SUNmat& wt)const;
  void maxTrSU3_2(SUNmat& gt,const SUNmat& wt)const;
  void maxTrSU3_3(SUNmat& gt,const SUNmat& wt)const;
#endif

  void max_trace_ovr(GaugeField1D& G,OverRelaxing OvR)const;
public:
  GaugeFixingStep(Condition cond,double or_prm,double sdm_prm,int Nvh)
    :Nvh_(Nvh),orp_(or_prm),sdmp_(sdm_prm){
    if(     cond == Coulomb) Ndim_= NDIM_-1;
    else if(cond == Landau)  Ndim_ = NDIM_;
    else abort();
  }

  const GaugeField1D umu_even(const GaugeField& Ue,const GaugeField& Uo)const;
  const GaugeField1D umu_odd( const GaugeField& Ue,const GaugeField& Uo)const;
  void step_SU2(GaugeField& Ue,GaugeField& Uo,OverRelaxing OvR)const;
  void step_sdm(GaugeField& Ue,GaugeField& Uo)const;

  void gauge_tr(GaugeField& U,const GaugeField1D& G)const;

  void gauge_tr_even(GaugeField& Ue,GaugeField& Uo,
		     const GaugeField1D& Ge)const;
  void gauge_tr_odd(GaugeField& Ue, GaugeField& Uo,
		    const GaugeField1D& Go)const;
};

#endif
