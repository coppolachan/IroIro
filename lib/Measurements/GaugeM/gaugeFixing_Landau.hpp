/*! @file gaugeFixing_Landau.hpp
    @brief definition of GaugeFixing_Landau 
    original version is written by H. Matsufuru
*/
#ifndef GAUGEFIXING_LANDAU_INCLUDED
#define GAUGEFIXING_LANDAU_INCLUDED

#include "gaugeFixing.hpp"
#include "gaugeFixingStep.hpp"
#include "include/pugi_interface.h"
#include "include/commonPrms.hpp"
#include "Geometry/shiftField.hpp"
#include <memory>
/*
  This class fix the gauge of configuration to the Landau gauge.
  Overrelaxation is incorporated.
  To escape the Gribov copy, if convergence is not reached within
  the iterations specified by Nreset, random gauge transformation
  is performed to reset the configuration.
  This is the reason that random number generator is needed at the
  construction of this class.
 */
class RandNum;

class GaugeFixing_Landau: public GaugeFixing{
 private:
  int    Niter_; // max iteration number   
  int    Nmeas_; // interval of measurements                         
  int    Nreset_;// Number of iteration to reset the config 
  int    Nor_;   // Number of iteration to start overrelaxation
  double esdm_;  // precision level to trigger the steepest descent method
  double prec_;  // convergence criterion  

  int Nvh_;

  const GaugeFixingStep* gstep_;
  const RandNum& rng_;

  double calc_F(const GaugeField& Ue,const GaugeField& Uo) const;
  double gauge_cond(const GaugeField& Ue,const GaugeField& Uo) const;
  void random_gtr(GaugeField& Ue,GaugeField& Uo)const;
  
 public:
  GaugeFixing_Landau(const RandNum& rng,XML::node GFnode)
    :rng_(rng),gstep_(NULL),Nvh_(CommonPrms::instance()->Nvol()/2){
    double orp,sdmp;
    //
    XML::read(GFnode,"max_iter",    Niter_, MANDATORY);
    XML::read(GFnode,"monitor_step",Nmeas_, MANDATORY);
    XML::read(GFnode,"reset_step",  Nreset_,MANDATORY);
    XML::read(GFnode,"or_step",     Nor_,   MANDATORY);
    XML::read(GFnode,"or_prm",      orp,    MANDATORY);
    XML::read(GFnode,"sdm_trigger", esdm_,  MANDATORY);
    XML::read(GFnode,"sdm_prm",     sdmp,   MANDATORY);
    XML::read(GFnode,"precision",   prec_,  MANDATORY);
    //
    gstep_= new GaugeFixingStep(Landau,orp,sdmp,Nvh_);
    if(!gstep_) abort();

    Mapping::init_shiftField_EvenOdd();
    CCIO::cout<<"GaugeFixing_Landau generated"<<std::endl;
  }

  GaugeFixing_Landau(const RandNum& rng,
		     int Niter,int Nmeas,int Nreset,
		     int Nor,double orp,double esdm,double sdmp,
		     double prec)
    :rng_(rng),Niter_(Niter),Nmeas_(Nmeas),Nreset_(Nreset),
     Nor_(Nor),esdm_(esdm),
     prec_(prec),Nvh_(CommonPrms::instance()->Nvol()/2),
     gstep_(new GaugeFixingStep(Landau,orp,sdmp,Nvh_)){
    //
    Mapping::init_shiftField_EvenOdd();
  }

  ~GaugeFixing_Landau(){delete gstep_;}

  const GaugeField do_fix(const GaugeField& Uin)const;
};

#endif  
