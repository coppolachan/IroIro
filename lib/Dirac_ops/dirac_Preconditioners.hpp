/*!
 * @file dirac_Preconditioners.hpp
 *
 * @brief Declaration of abstract class Preconditioner
 *
 * Concrete classes should be declared and specialized inside Dirac operators
 */
#include "include/field.h"

class Preconditioner {

public: 
  virtual const Field mult(const Field&) const = 0;    // pure virtual function
  virtual const Field mult_dag(const Field&) const = 0;// pure virtual function
  virtual const Field left(const Field&) const = 0;    // pure virtual function
  virtual const Field right(const Field&) const = 0;   // pure virtual function 
  virtual const Field left_dag(const Field&) const = 0;    // pure virtual function
  virtual const Field right_dag(const Field&) const = 0;   // pure virtual function 

};
