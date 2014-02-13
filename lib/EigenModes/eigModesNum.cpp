#include "eigModesNum.hpp"
#include "fopr_Linear.h"
#include "fopr_InvLinear.hpp"
#include "fopr_QuadLinear.h"
#include "fopr_RationalApprox.hpp"
#include "Solver/solver_CG.hpp"
#include "Solver/multiShiftSolver_CG.hpp"
#include "Tools/randNum_MP.h"
#include "Fields/field_expressions.hpp"

using namespace std;

void EigModesNum::count()const{
  count(cutoffs_);
}

void EigModesNum::count(const vector<double>& cutoff)const{
  using namespace FieldExpression;

  if(cutoff.empty()) abort();

  valarray<double> xi(DdagD_->fsize());

  CCIO::cout << setiosflags(ios_base::scientific);
  CCIO::cout<<"-----------------------------------------------------\n";
  CCIO::cout<<"       cutoff_squared      modes number below cutoff \n";
  CCIO::cout<<"-----------------------------------------------------\n";

  for(int i=0; i<cutoff.size(); ++i){

    double Lsq = cutoff[i]*cutoff[i];

    Fopr_Linear     DdDcut(1.0,Lsq,DdagD_);
    Solver_CG       slv(prec_,max_iter_,&DdDcut);
    Fopr_InvLinear  X1(-2.0*Lsq,1.0,&DdDcut,&slv);
    Fopr_QuadLinear X2(1.0,0.0,&X1);
    MultiShiftSolver_CG sslv(&X2,prec_,max_iter_);
    Fopr_RationalApprox Ra(r0_,res_,pole_,&X2,&sslv);
      
    double nu =0.0;
    for(int r=0; r<Nrand_; ++r){
      MPrand::mp_get_gauss(xi,*rng_);
      
      Field h2e(xi);
      Field XRa = X1.mult(Ra.mult(h2e));   
      
      h2e -= 2.0*XRa;
      h2e += X1.mult(Ra.mult(XRa)); /// (1 -X1*P(X2))*(1 -X1*P(X2))*xi
      nu += h2e*h2e;
    }
    nu /= 32.0*Nrand_;  /// 1/32 = 1/4^2*1/2, 1/4 is the prefactor of h2e
                        /// and 1/2 is the variant (= sigma^2) of xi
    
    CCIO::cout<<" "<<setw( 3)<<setiosflags(ios_base::right)<<i;
    CCIO::cout<<" "<<setw(10)<<setiosflags(ios_base::left) <<cutoff[i];
    CCIO::cout<<" "<<setw(20)<<setiosflags(ios_base::right)<<nu<<endl;
  }
}

void EigModesNum::init_RationalApprox(){
  RationalApprox_params RAprms(Ndeg_,Ndeg_,-1,2,100,epsilon_,1.0);
  RationalApprox RA(RAprms);
  r0_= RA.Const();
  res_= RA.Residuals();
  pole_= RA.Poles();
}
