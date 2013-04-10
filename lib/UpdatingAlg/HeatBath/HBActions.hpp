//--------------------------------------------------------------------
/*! @file HBActions.hpp
 *
 * @brief Definitions of Actions for Heat Bath updates
 *
 * Time-stamp: <2013-04-02 16:48:52 neo>
 */
//--------------------------------------------------------------------


/*!
  @file HBactions.hpp
  @brief Declaration of abstract HBAction class
 */
#ifndef HBACTION_INCLUDED
#define HBACTION_INCLUDED

#include <string>

#include "include/common_fields.hpp"

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
  virtual void hb_update(const RandNum&) = 0;

  ~HBAction(){}

};


#endif
