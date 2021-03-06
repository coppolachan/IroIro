/*!
  @file action.hpp

  @brief Declaration of abstract Action class, top level

  Time-stamp: <2013-04-16 11:23:55 neo>
 */

#ifndef ACTION_INCLUDED
#define ACTION_INCLUDED

#include <string>

#include "include/observer.hpp"
#include "include/common_fields.hpp"

class RandNum; // forward declaration 


/*!
 * @brief Definition of virtual Action class
 * Action class is defined 
 */
class Action : public Observer {
public:
  /*! Initializes the Action class */
  virtual void init(const RandNum&)= 0;

  /*!  Calculates action contribution */
  virtual double calc_H() = 0;

  /*! Calculates force contribution */
  virtual GaugeField md_force() = 0;

  /*! Default destructor */
  virtual ~Action(){}

  /*! Monitor force contribution */
  void monitor_force(GaugeField&, std::string);
};

enum {MetroStep, MDStep, PFStep};

#endif
