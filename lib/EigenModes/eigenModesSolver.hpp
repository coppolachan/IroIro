/*! @file eigenModesSolver.hpp
 *  @brief interface class calculate EigenModes 
 */
#ifndef EIGENMODESSOLVER_INCLUDED
#define EIGENMODESSOLVER_INCLUDED

#include <vector>
class Field;

class EigenModesSolver{
public:
  virtual void calc(std::vector<double>& lmd,std::vector<Field>& evec,int& N,
		    std::vector<Field>* exvec = NULL)const = 0;
};

#endif
