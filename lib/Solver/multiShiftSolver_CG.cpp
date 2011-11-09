/*!
 * @file multiShiftSolver_CG.cpp
 *
 * @brief Declaration of MultiShiftSolver_CG functions
 *
 */
#include "multiShiftSolver_CG.hpp"
#include <stdio.h>

using namespace std;

void MultiShiftSolver_CG::solve_init(vector<Field>& x,
				     vector<Field>& p,
				     Field& r,
				     Field& s,
				     double& rr,
				     vector<double>& zeta1,
				     vector<double>& zeta2,
				     vector<double>& csh2,
				     double& alphap,
				     double& betap) const{

  int Nshift = p.size();
  printf("number of shift = %d\n", Nshift);

  for(int i=0; i<Nshift; ++i){
    p[i] = s;
    x[i] = 0.0;
  }
  
  r = s;
  rr = r*r;  alphap = 0.0;
  betap  = 1.0;
}

void MultiShiftSolver_CG::solve_step(vector<Field>& x,
				     vector<Field>& p,
				     Field& r,
				     Field& s,
				     double& rr,
				     vector<double>& zeta1,
				     vector<double>& zeta2,
				     const vector<double>& sigma,
				     vector<double>& csh2,
				     double& alphap,
				     double& betap,
				     int& Nshift2, 
				     double& snorm, 
				     vector<double>& pp) const{

  using namespace FieldExpression;
  s = opr_->mult(p[0]);
  s += sigma[0]*p[0];

  double rrp = rr;
  double pap = s*p[0];
  double beta = -rrp/pap;

  x[0] -= beta*p[0];
  r += beta*s;
  rr = r*r;

  double alpha = rr/rrp;
  
  p[0] *= alpha;
  p[0] += r;
  pp[0] = rr;

  double alphah = 1.0 + alphap*beta/betap;
  for(int ish = 1; ish<Nshift2; ++ish){
    double zeta =(alphah-csh2[ish]*beta)/zeta1[ish]+(1.0-alphah)/zeta2[ish];
    zeta = 1.0/zeta;
    double zr = zeta/zeta1[ish];
    double betas  = beta  *  zr;
    double alphas = alpha * zr*zr;

    x[ish] -= betas * p[ish];
    p[ish] *= alphas;
    p[ish] += zeta * r;

    pp[ish] = p[ish] * p[ish];
    pp[ish] *= snorm;

    zeta2[ish] = zeta1[ish];
    zeta1[ish] = zeta;
  }

  for(int ish = Nshift2-1; ish>=0; --ish){
    //    printf("%4d %16.8e\n",ish,pp[ish]);
    if(pp[ish]> Params.GoalPrecision){
      Nshift2 = ish+1;
      break;
    }
  }
  alphap = alpha;
  betap  = beta;
}

void MultiShiftSolver_CG::solve(prop_t& xq, 
				const Field& b,
				const vector<double>& sigma, 
				double& diff,
				int& Nconv)const
{ 
  using namespace FieldExpression;
  
  cout << "Multi-shift solver Conjugate Gradient start" << endl;
  
  int Nshift = sigma.size();
  size_t fsize = b.size();
  
  printf("Number of shifts = %d\n", Nshift);
  printf(" -- values of sigma:\n");
  for(int i = 0; i<Nshift; ++i){ printf("   %d  %12.8f\n", i, sigma[i]);}
  
  double snorm = 1.0/b.norm();
  
  cout<<"fsize="<<b.size()<<" norm="<<b.norm()<<endl;
  Nconv = -1;
  
  Field s = b;
  Field r = b;
  vector<Field> p(Nshift);
  vector<Field> x(Nshift);

  vector<double> zeta1(Nshift,1.0);
  vector<double> zeta2(Nshift,1.0);
  vector<double> csh2(Nshift);
  vector<double> pp(Nshift);
  
  for(int i=0; i<Nshift; ++i){
    p[i].resize(fsize);
    x[i].resize(fsize);
    csh2[i] = sigma[i] -sigma[0];
  }
  
  double rr;
  double alphap, betap, rrp;
  int Nshift2 = Nshift;
  
  cout << "solve_init" << endl;
  solve_init(x,p,r,s,rr,zeta1,zeta2,csh2,alphap,betap);
  printf("  init: %22.15e\n",rr*snorm);
  
  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(x,p,r,s,rr,zeta1,zeta2,sigma,csh2,alphap,betap,
	       Nshift2,snorm,pp);
    printf("%6d  %22.15e  %4d\n",it,rr*snorm,Nshift2);
    
    if(rr*snorm < Params.GoalPrecision){
      Nconv = it;
      break;
    }
  }
  if(Nconv == -1) throw "Not converged.";
  
  printf("  residues of solutions:\n");
  diff = -1.0;
  for(int i=0; i<Nshift; ++i){
    s = opr_->mult(x[i]);
    s += sigma[i]*x[i];
    s -= b;
    double diff1 = s * s;
    diff1 *= snorm;
    printf("%6d  %22.15e\n",i,diff1);
    if(diff1>diff) diff = diff1;
  }
  printf("  diff(max) = %22.15e  \n",diff);
  
  for(int i=0; i<Nshift; ++i) xq[i] = x[i];
}
