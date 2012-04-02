/*! @file gaugeFixing_Landau.hpp
    @brief 
*/
#ifndef GAUGEFIXING_LANDAU_INCLUDED
#define GAUGEFIXING_LANDAU_INCLUDED

#include "include/pugi_interface.h"
#include "siteIndex_EvenOdd.hpp"
#include "include/common_field.h"
#include <valarray>
//! Landau gauge fixing.
/*
  This class fix the gauge of configuration to the Landau gauge.
  The algorithm is that developed by the Los Alamos group [see the
  implementation note].
  Overrelaxation is incorporated.
  To escape the Gribov copy, if convergence is not reached within
  the iterations specified by Nreset, random gauge transformation
  is performed to reset the configuration.
  This is the reason that random number generator is needed at the
  construction of this class.

  The implementation is not complete:
  - only applies to SU(3) case: because of specific implementation
      of maxTr function (Cabibbo-Marinari maximization).
  This should be improved in the version beyond test phase.
                                        [16 Feb 2012 H.Matsufuru]
 */
class RandNum;

class GaugeFixing_Landau{
 private:
  RandNum& rnd_;
  SiteIndex_EvenOdd* idx_;
  int    max_iter_;   // max iteration number
  int    Niter_;  // number of naive iterations
  int    Nmeas_;   // interval of measurements
  int    Nreset_;  // Number of iteration to reset the config.
  double prec_;   // convergence criterion
  double wp_;      // overrelaxation parameter

 public:
  GaugeFixing_Landau(const RandNum& rnd,XML::node GFnode):rnd_(rnd){
    XML::read(GFnode,"max_iteration",max_iter_);
    XML::read(GFnode,"naive_iteration",Niter_);
    XML::read(GFnode,"reset_iteration",Nreset_);
    XML::read(GFnode,"precision",prec_);
    XML::read(GFnode,"overrelax_param",wp_);
  }

  GaugeFixing_Landau(const RandNum& rnd,int Niter,int Nnaive,int Nmeas,
		     int Nreset,double enorm,double wp)
    :rnd_(rnd),idx_(SiteIndex_EvenOdd()::instance()),
     Niter_(Niter),Nmeas_(Nmeas),Nreset_(Nreset),
     enorm_(enorm),wp_(wp){}

  ~GaugeFixing_Landau(){}

  void fix(Field_G& Ufix, const Field_G& Uorg);

  void gauge_trf(Field_G& Ue, Field_G& Uo, Field_G& Geo, int Ieo);

  void set_randomGtrf(Field_G& Geo);

  //! one step of gauge fixing with overrelaxation parameter wp.
  void gfix_step(GaugeField& Ue, GaugeField& Uo, double wp);

  void calc_SG(double& sg,double& Fval,GaugeField& Ue,GaugeField& Uo);
  void calc_W(GaugeField& Weo,GaugeField& Ue,GaugeField& Uo, int Ieo);
  void calc_DLT(GaugeField& Weo,GaugeField& Ue,GaugeField& Uo, int Ieo);

  void maxTr(GaugeFieldField&,GaugeField&);
  void maxTr1(GaugeField&,GaugeField&);
  void maxTr2(GaugeField&,GaugeField&);
  void maxTr3(GaugeField&,GaugeField&);

};

#endif  
