/*!
  @file action.hpp
  
  @brief Declaration of abstract Action class

 */
#ifndef ACTION_INCLUDED
#define ACTION_INCLUDED

#include <string>

#include "include/macros.hpp"
#include "include/observer.hpp"


class Field;
class RandNum;

/*!
 * @brief Definition of virtual Action class
 * 
 * Action class is defined 
 */
class Action : public Observer {
public:
  /*!
    Initializes the Action class
   */
  virtual void init(const RandNum&, const void* x=0)= 0;

  /*!
    Calculates action contribution
   */
  virtual double calc_H()= 0;

  /*!
    Calculates force contribution
   */
  virtual Field md_force(const void* x=0)= 0;

  virtual ~Action(){}

  /*!
    Monitor force contribution
   */
  void monitor_force(Field&, std::string);

};

#endif
