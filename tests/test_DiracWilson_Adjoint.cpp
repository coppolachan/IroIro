#include "test_DiracWilson_Adjoint.hpp"
#include "Dirac_ops/dirac_wilson_adjoint.hpp"
#include "Tools/RandomNumGen/randNum_MT19937.h"
#include "Tools/randNum_MP.h"
#include "fopr.h"
#include "Tools/sunMatUtils.hpp"

using namespace std;
using namespace SUNmatUtils;

int Test_DiracWilson_Adjoint::run(){
  
  Field* u = &(input_.gconf->data);
  double m0 = 0.1;

  Dirac_Wilson_Adjoint Dwa(m0,u);

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  int length=4;
  RandNum_MT19937 rand(init, length);

  valarray<double>  ph(Dwa.fsize());
  MPrand::mp_get_gauss(ph,rand);
  
  Field psi(ph);
  Fopr_DdagD DdagD(&Dwa);
  Fopr_DDdag DDdag(&Dwa);
  //  Field DdDpsi = Dwa.mult_dag(Dwa.mult(psi));

  Field DdDpsi = DdagD.mult(psi);
  CCIO::cout<<" Re(psi*DdagD*psi)="<< psi*DdDpsi<<"\n";
  CCIO::cout<<" Im(psi*DdagD*psi)="<< psi.im_prod(DdDpsi)<<"\n";

  Field DDdpsi = DDdag.mult(psi);
  CCIO::cout<<" Re(psi*DDdag*psi)="<< psi*DDdpsi<<"\n";
  CCIO::cout<<" Im(psi*DDdag*psi)="<< psi.im_prod(DDdpsi)<<"\n";

  Field p5 = Dwa.gamma5(psi);
  CCIO::cout<<"psi.norm="<< psi.norm()<<" p5.norm()="<<p5.norm()<<"\n";

  Field Dpsi = Dwa.mult(psi);
  Field Ddpsi = Dwa.mult_dag(psi);
  CCIO::cout<<"Dpsi.norm="<< Dpsi.norm()<<" Ddpsi.norm()="<<Ddpsi.norm()<<"\n";

  for(int a=0;a<NADJ_;++a){
    for(int b=a+1;b<NADJ_;++b){
      
      SU3mat res = lambda_mul[a](lambda_mul[b](unity()));
      res -= lambda_mul[b](lambda_mul[a](unity()));
      
      SU3mat res1 
	= xI(lambda[su3alg[a*NADJ_+b]]()*2.0*su3str[a*NADJ_+b]);

      CCIO::cout<<"a="<<a<<" b="<<b<<"\n";      
      /*
      for(int i=0; i<NC_*NC_; ++i)
	CCIO::cout<<res.r(i) <<" "<<res.i(i) <<"     "
		  <<res1.r(i)<<" "<<res1.i(i)<<"\n";
      */
      for(int i=0; i<NC_*NC_; ++i)
	CCIO::cout<<res.r(i)-res1.r(i)<<" "<<res.i(i)-res1.i(i)<<"\n";
    }
  }

  return 0;
}
