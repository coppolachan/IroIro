/*!
 * @file solver_CG.cpp
 *
 * @brief Definition of Solver_CG class member functions
 *
 */
#include "Communicator/communicator.h"
#include "solver_CG.h"
using CommunicatorItems::pprintf;

using namespace std;

/*!
 * @brief Solves the linear equation
 *
 * @param xq  Output vector
 * @param b   Input vector
 * @param diff  Final residual
 * @param Nconv  Number of iterations for convergence
 */

void Solver_CG::solve(Field& xq, 
		      const Field& b, 
		      double& diff, 
		      int& Nconv) const{ 
  if(nodeid_==0) cout << "CG solver start" << endl;
  Nconv = -1;

  Field x = b;
  Field r = b;
  r -= opr_->mult(b);
  Field p = r;
  double rr = r*r;
  double snorm = b.norm();
  snorm = 1.0/snorm;

  if(Communicator::instance()->nodeid()==0)cout<<"snorm="<<snorm<<endl;
  if(nodeid_==0) pprintf("  init: %22.15e\n",rr*snorm);

  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(r,p,x,rr);
    if(nodeid_==0) pprintf("%6d  %22.15e\n",it,rr*snorm);
    
    if(rr*snorm < Params.GoalPrecision){
      Nconv = it;
      break;
    }
  }
  if(Nconv == -1) throw "Not converged.";

  p = opr_->mult(x);
  p -= b;
  diff = p.norm();

  xq = x;
}

inline void Solver_CG::solve_step(Field& r,Field& p,Field& x,double& rr)const {

  using namespace FieldExpression;

  Field s = opr_->mult(p);

  double pap = p*s;
  double rrp = rr;
  double cr = rrp/pap;

  Field v = p;
  x += cr*p;
  r -= cr*s;

  rr = r*r;
  p *= rr/rrp;
  p += r;
}
