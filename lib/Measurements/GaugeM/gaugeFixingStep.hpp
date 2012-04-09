/*! @file gaugeFixingStep.hpp
    @brief this class defines one step of the gauge fixing, which is 
    common to Coulomb and Landau gauge.
 */
#ifndef GAUGEFIXINGSTEP_INCLUDED
#define GAUGEFIXINGSTEP_INCLUDED

#include "include/common_fields.hpp"
#include "include/macros.hpp"

enum Condition{Coulomb,Landau};

class GaugeFixingStep{
private:
  const GaugeField1D max_trace(const GaugeField1D& W)const;
  void maxTrSU3_1(GaugeField1D& G,GaugeField1D& W)const;
  void maxTrSU3_2(GaugeField1D& G,GaugeField1D& W)const;
  void maxTrSU3_3(GaugeField1D& G,GaugeField1D& W)const;
  int Nvh_;
  int Ndim_;
public:
  GaugeFixingStep(Condition cndtn)
    :Nvh_(CommonPrms::instance()->Nvol()/2){
    if(     cndtn == Coulomb) Ndim_=NDIM_-1;
    else if(cndtn == Landau)  Ndim_ = NDIM_;
    else abort();
  }

  void gfix_step(GaugeField& Ue,GaugeField& Uo,double wp)const;
  
  void gauge_tr_even(GaugeField& Ue,GaugeField& Uo,const GaugeField1D& G)const;
  void gauge_tr_odd( GaugeField& Ue,GaugeField& Uo,const GaugeField1D& G)const;
  
  const GaugeField1D W_even(const GaugeField& Ue,const GaugeField& Uo)const;
  const GaugeField1D W_odd( const GaugeField& Ue,const GaugeField& Uo)const;
};

#endif
