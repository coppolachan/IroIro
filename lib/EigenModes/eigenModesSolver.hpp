/*! @file eigenModesSolver.hpp
 *  @brief interface class calculate EigenModes 
 */
#ifndef EIGENMODESSOLVER_INCLUDED
#define EIGENMODESSOLVER_INCLUDED

#include <vector>
class Field;

class EigenModesSolver{
public:
  virtual void calc(std::vector<double>&,std::vector<Field>&,int&)const = 0;
};

#endif
