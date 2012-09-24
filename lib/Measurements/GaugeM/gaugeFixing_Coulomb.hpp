/*! @file gaugeFixing_Coulomb.hpp
    @brief definition of GaugeFixing_Coulomb
    original version is written by H. Matsufuru
*/
#ifndef GAUGEFIXING_COULOMB_INCLUDED
#define GAUGEFIXING_COULOMB_INCLUDED

#include "gaugeFixing.hpp"
#include "gaugeFixingStep.hpp"
#include "include/pugi_interface.h"
#include "include/commonPrms.h"
#include "Main/Geometry/shiftField.hpp"
//#include "Main/Geometry/mapping.hpp"
#include "Tools/randNum_Factory.h"
#include <memory>
/*
  This class fix the gauge of configuration to the Coulomb gauge.
  Overrelaxation is incorporated.
  To escape the Gribov copy, if convergence is not reached within
  the iterations specified by Nreset, random gauge transformation
  is performed to reset the configuration.
  This is the reason that random number generator is needed at the
  construction of this class.
 */
class RandNum;

class GaugeFixing_Coulomb: public GaugeFixing{
 private:
  int    Niter_; // max iteration number
  int    Nmeas_; // interval of measurements
  int    Nreset_;// Number of iteration to reset the config
  int    Nor_;   // Number of iteration to start overrelaxation
  int    Nsdm_;  // Number of iteration to start steepest descent method
  double prec_;  // convergence criterion

  int Nvh_;

  const GaugeFixingStep* gstep_;
  //  std::auto_ptr<const RandNum> rng_;
  const RandNum& rng_;

  const std::vector<double> calc_F(const GaugeField& Ue,
				   const GaugeField& Uo)const;
  const std::vector<double> gauge_cond(const GaugeField& Ue,
				       const GaugeField& Uo)const;
  void random_gtr(GaugeField& Ue,GaugeField& Uo)const;

 public:
  GaugeFixing_Coulomb(const RandNum& rng, XML::node GFnode)
    :rng_(rng),gstep_(NULL),Nvh_(CommonPrms::instance()->Nvol()/2){
    double orp,sdmp;
    //
    XML::read(GFnode,"max_iter",    Niter_, MANDATORY);
    XML::read(GFnode,"monitor_step",Nmeas_, MANDATORY);
    XML::read(GFnode,"reset_step",  Nreset_,MANDATORY);
    XML::read(GFnode,"or_step",     Nor_,   MANDATORY);
    XML::read(GFnode,"or_prm",      orp,    MANDATORY);
    XML::read(GFnode,"sdm_step",    Nsdm_,  MANDATORY);
    XML::read(GFnode,"sdm_prm",     sdmp,   MANDATORY);
    XML::read(GFnode,"precision",   prec_,  MANDATORY);
    //
    gstep_= new GaugeFixingStep(Coulomb,orp,sdmp);
    if(!gstep_) abort();
    Mapping::init_shiftField_EvenOdd();
    CCIO::cout<<"GaugeFixing_Coulomb generated"<<std::endl;
  }

  GaugeFixing_Coulomb(const RandNum& rng, 
		      int Niter,int Nmeas,int Nreset,
		      int Nor,double orp,int Nsdm,double sdmp,
		      double prec)
    :rng_(rng),gstep_(new GaugeFixingStep(Coulomb,orp,sdmp)),
     Niter_(Niter),Nmeas_(Nmeas),Nreset_(Nreset),
     Nor_(Nor),Nsdm_(Nsdm),
     prec_(prec),Nvh_(CommonPrms::instance()->Nvol()/2){
    //
    Mapping::init_shiftField_EvenOdd();
  }

  ~GaugeFixing_Coulomb(){delete gstep_;}

  const GaugeField do_fix(const GaugeField& Uin)const;
};

#endif  
