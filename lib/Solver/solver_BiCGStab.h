/*!
 * @file solver_BiCGStab.h
 *
 * @brief Definition of BiCGstab solver
 *
 */
#ifndef SOLVER_BICGSTAB_INCLUDED
#define SOLVER_BICGSTAB_INCLUDED

#include <typeinfo>
#include "include/pugi_interface.h"
#include "include/fopr.h"
#include "Solver/solver.h"



struct Solver_BiCGStabParams{
  int MaxIter;
  double GoalPrecision;
  
  Solver_BiCGStabParams(XML::node Solver_node){
    XML::read(Solver_node, "MaxIter", MaxIter);
    XML::read(Solver_node, "Precision", GoalPrecision);
  }

  Solver_BiCGStabParams(const double prec_, const int max_){
    MaxIter = max_;
    GoalPrecision = prec_;
  }
};



class Solver_BiCGStab: public Solver {
private:
  const Fopr* opr_;
  const Solver_BiCGStabParams Params;
  mutable Field rh_;
  mutable Field s_;
  mutable Field t_;
  
  void solve_step(Field&, Field&, Field&, Field&, 
		  double&, double&, double&, double&) const;
public:
  Solver_BiCGStab(double precision, 
		  int MaxIterations,
		  const Fopr* fopr)
    :opr_(fopr),
     Params(Solver_BiCGStabParams(precision, MaxIterations)){
    rh_.resize(opr_->fsize());
    s_.resize(opr_->fsize());
    t_.resize(opr_->fsize());
  }
  
  Solver_BiCGStab(XML::node Solver_node,
		  const Fopr* fopr)
    :opr_(fopr),
     Params(Solver_BiCGStabParams(Solver_node)){
    rh_.resize(opr_->fsize());
    s_.resize(opr_->fsize());
    t_.resize(opr_->fsize());
  }
  
  
  ~Solver_BiCGStab(){}
  
  bool check_DdagD() const {
    return  (typeid(*opr_) == typeid(Fopr_DdagD));
  }

  void solve(Field& solution, 
	     const Field& source, 
	     double& diff, 
	     int& iterations) const;
};

#endif    
