/* 
 * @file solver.h
 *
 * @brief Definition of abstract Solver class and SolverOutput structure
 */
#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

#include <string>
#include "include/config.h" //defines VERBOSITY

class Field;

struct SolverOutput{
  double diff;   /*!< Residual */
  int Iterations; /*!< Iterations for convergence */
  std::string Msg;
  void print(std::string Msg2="");
};

class Solver{
public:
  virtual ~Solver(){}
  //virtual void set_mq(double) const = 0;
  //virtual void set_gconf(const Field&) const = 0;
  virtual SolverOutput solve(Field& solution, 
			     const Field& source) const =0;

  virtual bool check_DdagD() const  = 0;
};

#endif    
