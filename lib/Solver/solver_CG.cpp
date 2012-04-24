/*!
 * @file solver_CG.cpp
 *
 * @brief Definition of Solver_CG class member functions
 *
 */

#include <iostream>
#include <iomanip>
#include "Communicator/communicator.h"
#include "Fields/field_expressions.hpp"
#include "solver_CG.hpp"

using namespace std;

/*!
 * @brief Solves the linear equation (no preconditioning)
 *
 * @param xq  Output vector
 * @param b   Input vector
 */
SolverOutput Solver_CG::solve(Field& xq,const Field& b) const{ 

#if VERBOSITY>=SOLV_ITER_VERB_LEVEL
  CCIO::header("CG solver start");
#endif
  SolverOutput Out;
  Out.Msg = "CG solver";
  Out.Iterations = -1;

  TIMING_START;
  
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

  if(Out.Iterations == -1) {
    CCIO::cout<<" Not converged. Current residual: "<< rr*snorm << endl;
    abort();
  }

  p = opr_->mult(x);
  p -= b;
  Out.diff = p.norm();

  xq = x;

  TIMING_END(Out.timing);
    

  return Out;
}

inline void Solver_CG::solve_step(Field& r,Field& p,Field& x,double& rr)const {
  using namespace FieldExpression;
  double* x_ptr = x.getaddr(0);
  double* p_ptr = p.getaddr(0);
  double* r_ptr = r.getaddr(0);

  Field s = opr_->mult(p);//Ap

  double* s_ptr = s.getaddr(0);

  double pap = p*s;// (p,Ap)
  double rrp = rr;
  double cr = rrp/pap;// (r,r)/(p,Ap)

  // x = x + cr * p
  // r_k = r_k - cr * Ap
  // rr = (r_k,r_k)
  rr = 0.0;
  for (int i = 0; i < x.size(); ++i){
    x_ptr[i] += cr * p_ptr[i];
    r_ptr[i] -= cr * s_ptr[i];
    rr += r_ptr[i]*r_ptr[i];
  }
  rr = Communicator::instance()->reduce_sum(rr);
  cr = rr/rrp;
  for (int i = 0; i < p.size(); ++i){
    p_ptr[i] = p_ptr[i]*cr+r_ptr[i];
  }
  
  //p *= rr/rrp; // p = p*(r_k,r_k)/(r,r)
  //p += r; // p = p + p*(r_k,r_k)/(r,r)
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

  TIMING_START;
  
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

  TIMING_END(Out.timing);
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
