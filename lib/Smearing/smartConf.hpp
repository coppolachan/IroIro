/*!
  @file smartConf.hpp
  @brief Declares the SmartConf class
*/
#ifndef SMART_CONF_H_
#define SMART_CONF_H_

#include <complex>
#include "stoutSmear.hpp"
#include "include/observer.hpp"
#include "include/common_fields.hpp"

class Action;
/*!
  @brief Smeared configuration container

  It will behave like a configuration from the point of view of
  the HMC update and integrators.
  An "advanced configuration" object that can provide not only the 
  data to store the gauge configuration but also operations to manipulate
  it like smearing.
  
  It stores a list of smeared configurations.
*/
class SmartConf : public Observer {
  // Private members
  const int smearingLevels;
  Smear_Stout StoutSmearing;
  std::vector<GaugeField> SmearedSet;
  
  // Member functions
  void fill_smearedSet();
  GaugeField AnalyticSmearedForce(const GaugeField&, 
				  const GaugeField&) const;
  const GaugeField& get_smeared_conf(int) const;
  
  void set_iLambda(GaugeField& iLambda, 
		   GaugeField& e_iQ,
                   const GaugeField& iQ, 
		   const GaugeField& Sigmap,
                   const GaugeField& U)const;
  
  void set_uw(double& u, double& w,
              const SUNmat& iQ1, const SUNmat& iQ2)const ;
  void set_fj(std::complex<double>& f0, std::complex<double>& f1,
              std::complex<double>& f2, const double& u,
	      const double& w)const;
  
  double func_xi0(double w)const;
  double func_xi1(double w)const;
  
public:
  GaugeField* ThinLinks;      /*!< @brief Pointer to the thin 
				links configuration */

  /*! @brief XML constructor */
  //SmartConf(XML::node node, Field* Config, Smear_Stout& Stout)

  /*! @brief Standard constructor */
  SmartConf(int Nsmear, Smear_Stout& Stout):
    smearingLevels(Nsmear),
    StoutSmearing(Stout),
    ThinLinks(new GaugeField){
    for (int i=0; i< smearingLevels; ++i)
      SmearedSet.push_back(*(new GaugeField));
  }

  /*! For just thin links */
  SmartConf():smearingLevels(0),
	      StoutSmearing(),
	      SmearedSet(0),
	      ThinLinks(new GaugeField){}

  void set_GaugeField(){ fill_smearedSet(); }

  void smeared_force(GaugeField&) const;

  GaugeField* get_current_conf() const;

  GaugeField* select_conf(bool smeared) const {
    if (smeared){
      if (smearingLevels) return get_current_conf();
      else           	  return ThinLinks;
    }
    else return ThinLinks;
  }
  void observer_update() { fill_smearedSet();}
};

#endif
