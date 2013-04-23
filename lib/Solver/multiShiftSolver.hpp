/*! 
 * @file multiShiftSolver.hpp
 *
 * @brief Definition of abstract MultiShiftSolver class
 *
 * Time-stamp: <2013-04-23 11:30:39 neo>
 */

#ifndef MULTISOLVER_INCLUDED
#define MULTISOLVER_INCLUDED

#include "solver.hpp"

class Field;

class MultiShiftSolver{
public:
  virtual ~MultiShiftSolver(){}

  virtual SolverOutput solve(vector_Field& solution, 
			     const Field& source,
			     const vector_double& shifts,
			     double& diff,
			     int& Nconv) const = 0;
};

#endif   
