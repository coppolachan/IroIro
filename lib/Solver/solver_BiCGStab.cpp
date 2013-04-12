/*!
 * @file solver_BiCGStab.cpp
 *
 * @brief Definition of BiCGstab solver
 *
 */

#include <iostream>
#include <iomanip>

#include "solver_BiCGStab.hpp"
#include "Fields/field_expressions.hpp"
#include "Communicator/communicator.h"
#include "include/messages_macros.hpp"

using namespace std;

SolverOutput Solver_BiCGStab::
solve(Field& xq, const Field& b)const{ 
  using namespace FieldExpression;

#if VERBOSITY>=SOLV_ITER_VERB_LEVEL
  CCIO::header("BiCGStab solver start");
#endif
  SolverOutput Out;
  Out.Msg = "BiCGStab solver";
  Out.Iterations = -1;

  TIMING_START;

  double bnorm = b.norm();
  double snorm = 1.0/bnorm;

  _Message(SOLV_ITER_VERB_LEVEL, "b.norm()="<<bnorm<<std::endl);

  Field x = b;
  Field p(b.size());
  Field v(b.size());
  Field r = b;

  r-= opr_->mult(b);
  rh_= r;

  double rr = r*r;
  double rho = 1.0;
  double alp = 1.0;
  double omg = 1.0;

  _Message(SOLV_ITER_VERB_LEVEL," Init: "<< rr*snorm << std::endl);
  
  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(r,p,x,v,rr,rho,alp,omg);
#if VERBOSITY > SOLV_ITER_VERB_LEVEL
    CCIO::cout<< std::setw(5)<< "["<<it<<"] "
	      << std::setw(20) << rr*snorm<< "\n";
#endif
    if(rr*snorm < Params.GoalPrecision){
      Out.Iterations = it;
      break;
    }
  }
  if(Out.Iterations == -1) throw "Not converged.";

  p = opr_->mult(x) -b;
  Out.diff = p.norm();
  xq = x;

  TIMING_END(Out.timing);
  
  return Out;
}

void Solver_BiCGStab::
solve_step(Field& r, Field& p, Field& x, Field& v, 
	   double& rr, double& rho, double& alp, double& omg) const {

  using namespace FieldExpression;

  double rhop = rh_*r;
  double bet = rhop*alp/(rho*omg);
  
  p = r +bet*(p -omg*v);
  v = opr_->mult(p);
  alp = rhop/(rh_*v);

  s_= r -alp*v;
  t_= opr_->mult(s_);
  
  omg = (t_*s_)/(t_*t_);

  x += omg*s_ +alp*p;
  r = s_-omg*t_;

  rho = rhop;
  rr = r*r;
}
