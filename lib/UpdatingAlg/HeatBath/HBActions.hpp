//--------------------------------------------------------------------
/*! @file HBActions.hpp
 *
 * @brief Definitions of Actions for Heat Bath updates
 *
 */
//--------------------------------------------------------------------

//#include "Actions/action_gauge_wilson.hpp"

/*!
  @file action.hpp
  @brief Declaration of abstract Action class
 */
#ifndef HBACTION_INCLUDED
#define HBACTION_INCLUDED

#include <string>

#include "include/common_fields.hpp"
#include "Measurements/GaugeM/staples.hpp"
class RandNum;
/*!
 * @brief Definition of virtual HBAction class
 * Action class is defined for HeatBath
 */
class HBAction {
public:
  /*! Initializes the HB_Action class */
  virtual void init(const RandNum&)= 0;

  /*!  Calculates action contribution */
  virtual double calc_H() = 0;

  /*!  Updates action */
  virtual void hb_update() = 0;

  ~HBAction(){}

};


///////////////////////////////////////////////////


/*!
  Parameters for the ActionGaugeWilson class
 */
struct HBAction_SUNPrm {
  double beta;/*!< %Gauge coupling */

  HBAction_SUNPrm(const XML::node node){
    XML::read(node, "beta", beta);
  }
  HBAction_SUNPrm(const double beta_)
    :beta(beta_){}
};

/*!
  @brief Defines a concrete class for %gauge %Actions (Wilson)
  Specific for Heat-Bath updates
 */

class HBAction_SUN: public HBAction {
private:
  GaugeField* const u_;
  const Staples stpl_;
  HBAction_SUNPrm Params;
  const int Nvol_;


public:
  HBAction_SUN(const double beta, 
                    GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(HBAction_SUNPrm(beta)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  /*!
   * @brief Constructor with XML reader
   *
   */
  HBAction_SUN(const XML::node node,
                    GaugeField* const GField)
    :u_(GField),
     stpl_(),
     Params(HBAction_SUNPrm(node)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  ~HBAction_SUN(){}



  void init(const RandNum&){};
  double calc_H(){};
  void hb_update(){};

  
};

#endif
