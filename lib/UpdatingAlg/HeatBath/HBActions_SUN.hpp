//--------------------------------------------------------------------
/*! @file HBActions.hpp
 *
 * @brief Definitions of actions for Heat Bath updates
 *
 * Time-stamp: <2013-04-03 14:33:30 neo>
 */
//--------------------------------------------------------------------
#ifndef HBACTION_SUN_INCLUDED
#define HBACTION_SUN_INCLUDED


#include "HBActions.hpp"
#include "Measurements/GaugeM/staples.hpp"

///////////////////////////////////////////////////


/*!
  Parameters for the ActionGaugeWilson class
 */
struct HBAction_SUNPrm {
  double beta;/*!< %Gauge coupling */
  double beta_adj;/*!< %Gauge coupling */


  HBAction_SUNPrm(const XML::node node){
    XML::read(node, "beta", beta);
    XML::read(node, "beta_adj", beta_adj);
  }
  HBAction_SUNPrm(const double beta_, const double beta_adj_)
    :beta(beta_),
     beta_adj(beta_adj_){}
};



/*!
  @brief Defines a concrete class for %gauge %Actions (Wilson)
  Specific for Heat-Bath updates
 */

class HBAction_SUN: public HBAction {
private:
  GaugeField* const u_;
  const Staples stplfield_;
  HBAction_SUNPrm Params_;
  const int Nvol_;


public:
  HBAction_SUN(const double beta, const double beta_adj,
                    GaugeField* const GField)
    :u_(GField),
     stplfield_(),
     Params_(HBAction_SUNPrm(beta, beta_adj)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  /*!
   * @brief Constructor with XML reader
   *
   */
  HBAction_SUN(const XML::node node,
                    GaugeField* const GField)
    :u_(GField),
     stplfield_(),
     Params_(HBAction_SUNPrm(node)), 
     Nvol_(CommonPrms::instance()->Nvol()){}
  
  ~HBAction_SUN(){}



  void init(const RandNum&);
  double calc_H();
  void hb_update(const RandNum&);

  
};


#endif
