/*!
 * @file solver_CG.cpp
 *
 * @brief Definition of Solver_CG class member functions
 *
 */

#include <iostream>
#include <iomanip>
#include "Communicator/communicator.h"
#include "solver_CG.h"
using CommunicatorItems::pprintf;

using namespace std;

/*!
 * @brief Solves the linear equation (no preconditioning)
 *
 * @param xq  Output vector
 * @param b   Input vector
 */
SolverOutput Solver_CG::solve(Field& xq, 
		      const Field& b) const{ 
#if VERBOSITY>1
  CCIO::header("CG solver start");
#endif
  SolverOutput Out;
  Out.Msg = "CG solver";
  Out.Iterations = -1;
  
  Field x = b;//initial condition
  Field r = b;//initial residual
  r -= opr_->mult(x);
  Field p = r;
  double rr = r*r;// (r,r)
  double snorm = b.norm();
  snorm = 1.0/snorm;

#if VERBOSITY>1
  CCIO::cout<<" Snorm = "<< snorm << endl;
  CCIO::cout<<" Init  = "<< rr*snorm<< endl;
#endif
  for(int it = 0; it < Params.MaxIter; ++it){
    solve_step(r,p,x,rr);
#if VERBOSITY>1
    CCIO::cout<< std::setw(5)<< "["<<it<<"] "
	      << std::setw(20) << rr*snorm<< "\n";
#endif    
    if(rr*snorm < Params.GoalPrecision){
      Out.Iterations = it;
      break;
    }
  }
  if(Out.Iterations == -1) throw "Not converged.";

  p = opr_->mult(x);
  p -= b;
  Out.diff = p.norm();

  xq = x;

  return Out;
}

inline void Solver_CG::solve_step(Field& r,Field& p,Field& x,double& rr)const {

  using namespace FieldExpression;

  Field s = opr_->mult(p);//Ap

  double pap = p*s;// (p,Ap)
  double rrp = rr;
  double cr = rrp/pap;// (r,r)/(p,Ap)

  x += cr*p; // x = x + cr * p
  r -= cr*s; // r_k = r_k - cr * Ap

  rr = r*r; // rr = (r_k,r_k)
  p *= rr/rrp; // p = p*(r_k,r_k)/(r,r)
  p += r; // p = p + p*(r_k,r_k)/(r,r)
}

//////////////////////////////////////////////////////////////////////////////////
/*!
 * @brief Solves the linear equation (preconditioned version)
 *
 * @param xq  Output vector
 * @param b   Input vector
 */
SolverOutput Solver_CG_Precondition::solve(Field& xq, 
				   const Field& b) const{ 

#if VERBOSITY > 1
  CCIO::header("CG solver Preconditioned start");
#endif
  SolverOutput Out;
  Out.Msg = "CG solver Preconditioned";
  Out.Iterations = -1;
  
  Field x = b;//starting guess
  Field r = b;
  r -= opr_->mult_prec(x);
  Field p = r;
  double rr = r*r;
  double snorm = b.norm();
  snorm = 1.0/snorm;

#if VERBOSITY > 1
  CCIO::cout<<" Snorm = "<<snorm<<std::endl;
  CCIO::cout<<" Init  = "<< rr*snorm << std::endl;
#endif

  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(r,p,x,rr);
#if VERBOSITY > 1
    CCIO::cout<< std::setw(5)<< "["<<it<<"] "<< std::setw(20) << rr*snorm<< "\n";
#endif    
  
    if(rr*snorm < Params.GoalPrecision){
      Out.Iterations = it;
      break;
    }
  }
  if(Out.Iterations == -1) throw "Not converged.";

  p = opr_->mult_prec(x);
  p -= b;
  Out.diff = p.norm();

  xq = x;
  return Out;
}

inline void Solver_CG_Precondition::solve_step(Field& r,Field& p,Field& x,double& rr)const {
  
  using namespace FieldExpression;
  
  Field s = opr_->mult_prec(p);
  
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
