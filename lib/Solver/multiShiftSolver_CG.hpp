/*!
 * @file multiShiftSolver_CG.h
 *
 * @brief Definition of the MultiShiftSolver_CG class
 *
 */

#ifndef MULTISHIFTSOLVER_CG_INCLUDED
#define MULTISHIFTSOLVER_CG_INCLUDED

#include "include/pugi_interface.h"
#include "include/fopr.h"
#include "multiShiftSolver.h"

typedef std::vector<Field> prop_t;

struct MultiShiftSolver_CG_Params{
  double GoalPrecision;
  int MaxIter;

  MultiShiftSolver_CG_Params(XML::node MShift_node){
    XML::read(MShift_node, "MaxIter", MaxIter);
    XML::read(MShift_node, "Precision", GoalPrecision);
  }

  MultiShiftSolver_CG_Params(double prec_, int iter_){
    GoalPrecision = prec_;
    MaxIter = iter_;
  }
};

class MultiShiftSolver_CG : public MultiShiftSolver {

private:
  const Fopr* opr_;
  const MultiShiftSolver_CG_Params Params;

  void solve_init(std::vector<Field>&,
		  std::vector<Field>&,
		  Field&,
		  Field&, 
		  double&,
		  std::vector<double>&, 
		  std::vector<double>&,
		  std::vector<double>&,
		  double&, 
		  double&) const;

  void solve_step(std::vector<Field>&, 
		  std::vector<Field>&,
		  Field&, 
		  Field&, 
		  double&,
		  std::vector<double>&,
		  std::vector<double>&,
		  const std::vector<double>&,
		  std::vector<double>&,
		  double&,
		  double&, 
		  int&, 
		  double&, 
		  std::vector<double>&) const;

public:
  MultiShiftSolver_CG(const Fopr* fopr, 
		      const double stop_cond,
		      const int Niter)
    :opr_(fopr),
     Params(MultiShiftSolver_CG_Params(stop_cond, Niter)){}

  MultiShiftSolver_CG(const Fopr* fopr, 
		      const XML::node node)
    :opr_(fopr),
     Params(MultiShiftSolver_CG_Params(node)){}

  ~MultiShiftSolver_CG(){}
  
  void solve(prop_t& solution, 
	     const Field& source,
	     const std::vector<double>& shifts,
	     double& residual,
	     int& Nconv) const;

};

#endif    
