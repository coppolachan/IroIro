/*!
 * @file multiShiftSolver_CG.hpp
 * @brief Definition of the MultiShiftSolver_CG class
 * Time-stamp: <2014-07-02 21:33:57 noaki>
 */
#ifndef MULTISHIFTSOLVER_CG_INCLUDED
#define MULTISHIFTSOLVER_CG_INCLUDED

#include "include/pugi_interface.h"
#include "Fopr/fopr.h"
#include "multiShiftSolver.hpp"

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

  void solve_init(vector_Field&,
		  vector_Field&,
		  Field&,
		  Field&, 
		  double&,
		  vector_double&, 
		  vector_double&,
		  vector_double&,
		  double&, 
		  double&) const;

  void solve_step(vector_Field&, 
		  vector_Field&,
		  Field&, 
		  Field&, 
		  double&,
		  vector_double&,
		  vector_double&,
		  const vector_double&,
		  vector_double&,
		  double&,
		  double&, 
		  int&, 
		  double&, 
		  vector_double&) const;

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
  
  SolverOutput solve(vector_Field& solution, 
		     const Field& source,
		     const vector_double& shifts,
		     double& residual,
		     int& Nconv) const;

};

#endif    
