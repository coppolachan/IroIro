#include "eigModesNum.hpp"
#include "Fopr/fopr_Linear.h"
#include "Fopr/fopr_InvLinear.hpp"
#include "Fopr/fopr_QuadLinear.h"
#include "Fopr/fopr_Chebyshev.h"
#include "Solver/solver_CG.hpp"
#include "Tools/randNum_MP.h"
#include "Fields/field_expressions.hpp"
#include <gsl/gsl_integration.h>

using namespace std;

void EigModesNum::do_count()const{
  if(cutoffs_.empty()) abort();
  CCIO::cout<<"stochastic eigenmode-counting with "
	    << Nrand_<<" sources and "<<cutoffs_.size()<<" lambda's\n";

  CCIO::cout << setiosflags(ios_base::scientific);
  CCIO::cout<<"----------------------------------------------------\n";
  CCIO::cout<<"      cutoff (lambda)     modes number below cutoff \n";
  CCIO::cout<<"----------------------------------------------------\n";

  vector<double> nu(cutoffs_.size(),0.0);
  valarray<double> xi(DdagD_->fsize());

  for(int r=0; r<Nrand_; ++r){
    MPrand::mp_get_gauss(xi,*rng_);
    CCIO::cout<<"noise source "<<r<<"\n";
    
    for(int i=0; i<cutoffs_.size(); ++i){
      Field f(xi);
      hproj(f,cutoffs_[i]*cutoffs_[i]/Mratio_/Mratio_);
      double nu_r = f*f;
      nu[i] += nu_r;

      CCIO::cout<<" "<<setw( 3)<<setiosflags(ios_base::right)<< i;
      CCIO::cout<<" "<<setw(10)<<setiosflags(ios_base::left )<< cutoffs_[i];
      CCIO::cout<<" "<<setw(20)<<setiosflags(ios_base::right)<< nu_r <<endl;
    }
  }
  CCIO::cout<<"average over sources:\n";
  for(int i=0; i<cutoffs_.size(); ++i){
    CCIO::cout<<" "<<setw( 3)<<setiosflags(ios_base::right)<< i;
    CCIO::cout<<" "<<setw(10)<<setiosflags(ios_base::left )<< cutoffs_[i];
    CCIO::cout<<" "<<setw(20)<<setiosflags(ios_base::right)<< nu[i]/Nrand_ <<endl;
  }
}

void EigModesNum::hproj(Field& f,double Msq)const{
  using namespace FieldExpression;

  Fopr_Linear     DdDcut(1.0,Msq,DdagD_);
  Solver_CG       slv(prec_,max_iter_,&DdDcut);
  Fopr_InvLinear  X1(-2.0*Msq,1.0,&DdDcut,&slv);
    
  Fopr_QuadLinear X2(2.0/(1.0-epsilon_),
		     (-1.0-epsilon_)/(1.0-epsilon_),&X1);
  Fopr_Chebyshev  PX2(&X2,coeff_);

  //// evaluation of (1 -X1*P(X2))*(1 -X1*P(X2))*f
  Field xPx2 = X1.mult(PX2.mult(f)); /// X1*P(X2)*f
  f -= 2.0*xPx2;                     
  f += X1.mult(PX2.mult(xPx2));    /// (X1*P(X2))**2*f
  f /= 4.0;
  
  /* /// check of the functional form
  int N=1000;
  for(int i=-N; i<N; ++i){
    double y = i/double(N);

    double sgn = func_sgn(-2.0*Msq/(y*y +Msq)+1.0,epsilon_,coeff_);
    double sgn_fopr = X1.func(y)*PX2.func(y);
    CCIO::cout<<y<<" "<<sgn<<"  "<<sgn_fopr<<"\n";
  }
  */
}

double EigModesNum::func_sgn(double y,double eps,const vector<double>& c){
  return y*Chebyshev::series((2.0*y*y-1.0-eps)/(1.0-eps),c);  
}

double EigModesNum::func_rh4(double y,void* prms){
  vector<double> coeff(*static_cast<vector<double>*>(prms));
  double eps = coeff.back();
  coeff.pop_back();
  double hx = 0.5*(1.0-func_sgn(y,eps,coeff));
  return (1.0 +y)/sqrt(1.0-y*y)/(1.0-y*y)*hx*hx*hx*hx;
}

double EigModesNum::func_sqrtInv(double x,void* prms){
  return 1.0/sqrt(x);}

void EigModesNum::init_ChebyshevApprox(){

  /// Below, "100" is the order of the Chebyshev approx. 
  /// We truncate the resulting series to Npoly_'th order(Npoly_<=100).
  Chebyshev::approx(coeff_,100,func_sqrtInv,epsilon_,1.0);

  CCIO::cout<<"coefficients of the Chebyshev series:\n";
  for(int k=0; k<coeff_.size(); ++k)
    CCIO::cout<<k<<" "<<coeff_[k]<<"\n";  
  CCIO::cout<<"-------------------------------------\n";

  //// estimation of delta
  int N= 1000;
  double delta = 0.0;
  for(int i=0; i<=N; ++i){
    double y = epsilon_+i/double(N)*(1.0-epsilon_);
    double z = (2.0*y-1.0-epsilon_)/(1.0-epsilon_);
    double diff = fabs(1.0-sqrt(y)*Chebyshev::series(z,coeff_));
    //CCIO::cout<<"y="<<y<<" diff= "<<diff<<"\n";
    delta = max(delta,diff);
  }
  CCIO::cout<<"delta= "<<delta<<"\n";

  //// estimation of M/Mster 
  double res, err;
  gsl_function F;  
  vector<double> param(coeff_);
  param.push_back(epsilon_);

  F.function = &EigModesNum::func_rh4;
  F.params = &param;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(200);
  gsl_integration_qag(&F,-sqrt(epsilon_),sqrt(epsilon_),
		      0.0,1e-7,200,2,w,&res,&err);
  
  Mratio_= res+ sqrt((1.0-sqrt(epsilon_))/(1.0+sqrt(epsilon_)));

  CCIO::cout<<"M/Mstar="<<Mratio_<<"\n";
  gsl_integration_workspace_free(w);

  /*
  //// actual curves of sign function with the current approximation
  for(int i=-N; i<N; ++i){
    double y = i/double(N);
    double sgn = func_sgn(y,epsilon_,coeff_);
    CCIO::cout<<y<<" "<<sgn<<"\n";
  }
  */
}
