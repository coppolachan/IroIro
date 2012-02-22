/* 
 * @file multiShiftSolver.h
 *
 * @brief Definition of abstract MultiShiftSolver class
 */
#ifndef MULTISOLVER_INCLUDED
#define MULTISOLVER_INCLUDED


#include <vector>

#include "Fields/field_expressions.hpp"

class Field;

class MultiShiftSolver{
public:
  virtual ~MultiShiftSolver(){}
  //virtual void set_mq(double) const = 0;
  //virtual void set_gconf(const Field&) const = 0;
  virtual void solve(std::vector<Field>& solution, 
		     const Field& source,
		     const std::vector<double>& shifts,
		     double& diff,
		     int& Nconv) const = 0;
};

#endif   
