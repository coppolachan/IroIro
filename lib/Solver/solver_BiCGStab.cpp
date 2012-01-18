/*!
 * @file solver_BiCGStab.cpp
 *
 * @brief Definition of BiCGstab solver
 *
 */
#include "solver_BiCGStab.h"
#include "Fields/field_expressions.hpp"
#include "Communicator/communicator.h"

using CommunicatorItems::pprintf;
using namespace std;

SolverOutput Solver_BiCGStab::
solve(Field& xq, const Field& b)const{ 
  using namespace FieldExpression;

#if VERBOSITY > 1
  CCIO::header("BiCGStab solver start");
#endif
  SolverOutput Out;
  Out.Msg = "BiCGStab solver";
  Out.Iterations = -1;


  double bnorm = b.norm();
  double snorm = 1.0/bnorm;

#if VERBOSITY > 1
  CCIO::cout<<"b.norm()="<<bnorm<<std::endl;
#endif

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

#if VERBOSITY > 1
  CCIO::cout << " Init: "<< rr*snorm << std::endl;
#endif
  
  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(r,p,x,v,rr,rho,alp,omg);
#if VERBOSITY > 1
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
