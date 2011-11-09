/* 
 * @file solver.h
 *
 * @brief Definition of abstract Solver class
 */
#ifndef SOLVER_INCLUDED
#define SOLVER_INCLUDED

class Field;

class Solver{
public:
  virtual ~Solver(){}
  //virtual void set_mq(double) const = 0;
  //virtual void set_gconf(const Field&) const = 0;
  virtual void solve(Field& solution, 
		     const Field& source,
		     double& diff,
		     int& nconv) const =0;

  virtual bool check_DdagD() const  = 0;
};

#endif    
