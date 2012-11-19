/*!
 * @file multiShiftSolver_CG.cpp
 *
 * @brief Declaration of MultiShiftSolver_CG functions
 *
 */
#include "multiShiftSolver_CG.hpp"
#include "Fields/field_expressions.hpp"
#include "include/messages_macros.hpp"

typedef struct FermionSpinor{
  double _Complex v[12];
}Spinor;

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
  _Message(SOLV_ITER_VERB_LEVEL, "    MultiShiftSolver_CG Inizialitation\n");

  for(int i=0; i<Nshift; ++i){
    p[i] = s;
    x[i] = 0.0;
  }
  
  r = s;
  rr = r*r;  alphap = 0.0;
  betap  = 1.0;
  _Message(SOLV_ITER_VERB_LEVEL, "    | Initial residual |^2 = "<<rr<<"\n");
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
    if(pp[ish]> Params.GoalPrecision){
      Nshift2 = ish+1;
      break;
    }
  }
  alphap = alpha;
  betap  = beta;
}

SolverOutput MultiShiftSolver_CG::solve(prop_t& xq, 
				const Field& b,
				const vector<double>& sigma, 
				double& diff,
				int& Nconv)const
{ 
  using namespace FieldExpression;
  
  _Message(SOLV_ITER_VERB_LEVEL, "Multi-shift solver Conjugate Gradient start\n");

  SolverOutput Out;
  Out.Msg = "Multishift CG solver";
  Out.Iterations = -1;

  TIMING_START;

  int Nshift = sigma.size();
  size_t fsize = b.size();

  double snorm = 1.0/b.norm();
  Nconv = -1;
  
  Field s = b;
  Field r = b;
  vector<Field> p(Nshift);
  vector<Field> x(Nshift);

  vector<double> zeta1(Nshift,1.0);
  vector<double> zeta2(Nshift,1.0);
  vector<double> csh2(Nshift);
  vector<double> pp(Nshift);

  double rr;
  double alphap, betap, rrp;
  int Nshift2 = Nshift;

  // Initial messages
  _Message(SOLV_ITER_VERB_LEVEL, "    -------------------\n");
  _Message(SOLV_ITER_VERB_LEVEL, "    Number of shifts = "<< Nshift<<"\n");
  _Message(SOLV_ITER_VERB_LEVEL, "    Values of shifts:\n");
  for(int i = 0; i<Nshift; ++i){
    _Message(SOLV_ITER_VERB_LEVEL, "      #["<<i<<"] = "<< sigma[i]<<"\n");
  }
  _Message(SOLV_ITER_VERB_LEVEL, "    -------------------\n");

  // Initial condition
  for(int i=0; i<Nshift; ++i){
    p[i].resize(fsize);
    x[i].resize(fsize);
    csh2[i] = sigma[i] -sigma[0];
  }
  
  solve_init(x,p,r,s,rr,zeta1,zeta2,csh2,alphap,betap);
  _Message(SOLV_ITER_VERB_LEVEL, "    | Init | = "<<rr*snorm<<"\n");
  
  for(int it = 0; it < Params.MaxIter; it++){
    solve_step(x,p,r,s,rr,zeta1,zeta2,sigma,csh2,alphap,betap,
	       Nshift2,snorm,pp);
    double residual = rr*snorm; 
    _Message(SOLV_ITER_VERB_LEVEL, "   "<<std::setw(5)<<"["<<it<<"]  "
	     <<std::setw(20)<<residual<<"     Left: "<<Nshift2<<" \n");
    
    if(residual < Params.GoalPrecision){
      Nconv = it;
      break;
    }
  }
  if(Nconv == -1) throw "Not converged.";
  
  _Message(SOLV_ITER_VERB_LEVEL, "  --- Summary of true residuals\n");
  diff = -1.0;
  for(int i=0; i<Nshift; ++i){
    s = opr_->mult(x[i]);
    s += sigma[i]*x[i];
    s -= b;
    double diff1 = s * s;
    diff1 *= snorm;
    _Message(SOLV_ITER_VERB_LEVEL, "       ["<<i<<"]  "<<diff1<<"\n");
    if(diff1>diff) diff = diff1;
  }
  _Message(SOLV_ITER_VERB_LEVEL, " Maximum residual  = "<<diff<<"\n");
  
  for(int i=0; i<Nshift; ++i) xq[i] = x[i];

  Out.diff = sqrt(diff);
  TIMING_END(Out.timing);
  
 return Out;
}
