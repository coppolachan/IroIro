/*! @file gaugeFixing_Landau.hpp
    @brief definition of GaugeFixing_Landau 
    original version is written by H. Matsufuru
*/
#ifndef GAUGEFIXING_LANDAU_INCLUDED
#define GAUGEFIXING_LANDAU_INCLUDED

#include "gaugeFixing.hpp"
#include "gaugeFixingStep.hpp"
#include "include/pugi_interface.h"
#include "include/commonPrms.h"
#include "Main/Geometry/mapping.hpp"
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
  int    max_iter_; // max iteration number
  int    Niter_;    // number of naive iterations
  int    Nmeas_;    // interval of measurements
  int    Nreset_;   // Number of iteration to reset the config.
  double prec_;     // convergence criterion
  double wp_;       // overrelaxation parameter
  int Nvh_;

  const RandNum& rnd_;
  const GaugeFixingStep* gstep_;

  double calc_F(const GaugeField& Ue,const GaugeField& Uo) const;
  double calc_SG(const GaugeField& Ue,const GaugeField& Uo) const;
  void random_gtr(GaugeField& Ue,GaugeField& Uo)const;
  
 public:
  GaugeFixing_Landau(const RandNum& rnd,XML::node GFnode)
    :rnd_(rnd),gstep_(new GaugeFixingStep(Landau)){
    XML::read(GFnode,"max_iteration",max_iter_);
    XML::read(GFnode,"naive_iteration",Niter_);
    XML::read(GFnode,"reset_iteration",Nreset_);
    XML::read(GFnode,"precision",prec_);
    XML::read(GFnode,"overrelax_param",wp_);
    //
    Mapping::init_shiftField_EvenOdd();
  }

  GaugeFixing_Landau(const RandNum& rnd,int Niter,int Nnaive,int Nmeas,
		     int Nreset,double prec,double wp)
    :rnd_(rnd),gstep_(new GaugeFixingStep(Landau)),
     Niter_(Niter),Nmeas_(Nmeas),Nreset_(Nreset),prec_(prec),wp_(wp),
     Nvh_(CommonPrms::instance()->Nvol()/2){
    //
    Mapping::init_shiftField_EvenOdd();
  }

  ~GaugeFixing_Landau(){delete gstep_;}

  const GaugeField fix(const GaugeField& Uin)const;
};

#endif  
