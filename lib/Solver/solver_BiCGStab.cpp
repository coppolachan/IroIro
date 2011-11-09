/*!
 * @file solver_BiCGStab.cpp
 *
 * @brief Definition of BiCGstab solver
 *
 */
#include "solver_BiCGStab.h"
#include "Communicator/communicator.h"

using CommunicatorItems::pprintf;
using namespace std;

void Solver_BiCGStab::
solve(Field& xq, const Field& b, double& diff, int& Nconv)const{ 
  using namespace FieldExpression;

  if(Communicator::instance()->nodeid()==0)
    cout<<"BiCGStab solver start"<<endl;

  double bnorm = b.norm();
  double snorm = 1.0/bnorm;
  if(Communicator::instance()->nodeid()==0) cout<<"b.norm()="<<bnorm<<endl;
  Nconv = -1;

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

  pprintf("solve_init");
  pprintf("  init: %22.15e\n",rr*snorm);

  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(r,p,x,v,rr,rho,alp,omg);
    pprintf("%6d  %22.15e\n",it,rr*snorm);

    if(rr*snorm < Params.GoalPrecision){
      Nconv = it;
      break;
    }
  }
  if(Nconv == -1) throw "Not converged.";

  p = opr_->mult(x) -b;
  diff = p.norm();
  xq = x;
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
