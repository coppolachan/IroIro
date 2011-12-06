/*!
  @file action.h
  
  @brief Declaration of abstract Action class

 */
#ifndef ACTION_INCLUDED
#define ACTION_INCLUDED

class Field;
class RandNum;

/*!
 * @class Action class
 * @brief Definition of virtual Action class
 * 
 * Action class is defined 
 */
class Action{
public:
  /*!
    Initializes the Action class
   */
  virtual void init(Field&,const RandNum&, const void* x=0)= 0;

  /*!
    Calculates action contribution
   */
  virtual double calc_H()= 0;

  virtual Field md_force(const void* x=0)= 0;

  virtual ~Action(){};
};

#endif
