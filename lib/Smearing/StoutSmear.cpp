/*
  @file StoutSmear.hpp

  @brief Defines Stout smearing class

*/

#include "StoutSmear.hpp"

using namespace std;

#include "Communicator/comm_io.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/messages_macros.hpp"
#include <complex>


typedef complex<double> dcomplex;

//====================================================================
void Smear_Stout::smear(GaugeField& u_smr, const GaugeField& u_in) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Nvol = CommonPrms::Nvol();

  GaugeField u_tmp1, q_mu;
  GaugeField1D U_mu;
  SUNmat ut;

  _Message(DEBUG_VERB_LEVEL, "Stout smearing started\n");

  SmearBase->smear(u_tmp1,u_in);

  for(int mu = 0; mu < NDIM_; ++mu){
    U_mu = DirSlice(u_in, mu);
    for(int site = 0; site < Nvol; ++site){
      ut = matrix(u_tmp1,site,mu) * matrix_dag(U_mu,site);//Omega_mu
      SetMatrix(q_mu, anti_hermite_traceless(ut), site, mu);
    }
  }



  exponentiate_iQ(u_tmp1, q_mu);
  
  for(int mu = 0; mu < NDIM_; ++mu){
    U_mu = DirSlice(u_in, mu);
    for(int site = 0; site < Nvol; ++site){
      ut = matrix(u_tmp1,site,mu) * matrix(U_mu,site);//New smeared matrix
      SetMatrix(u_smr, ut, site, mu);
    }
  }

  _Message(DEBUG_VERB_LEVEL, "Stout smearing completed \n");
  
}
//====================================================================
void Smear_Stout::BaseSmear(GaugeField& C, const GaugeField& u_in) const{
  SmearBase->smear(C, u_in);
}
//====================================================================
void Smear_Stout::derivative(GaugeField& SigmaTerm, 
			     const GaugeField& iLambda,
			     const GaugeField& Gauge) const{
  
  SmearBase->derivative(SigmaTerm, iLambda, Gauge);
}
//====================================================================
void Smear_Stout::exponentiate_iQ(GaugeField& e_iQ, const GaugeField& iQ) const{
  using namespace SUNmatUtils;
  using namespace FieldUtils;

  int Nvol = e_iQ.format.Nvol();
  int Nex  = e_iQ.format.Nex();

  SUNmat iQ0, iQ1, iQ2, iQ3;
  iQ0 = unity();

  double c0, c1, c0max, u_val, w, theta, xi0, u2, w2, cosw, fden;
  dcomplex f0, f1, f2, h0, h1, h2, e2iu, emiu, ixi0, qt;

  for(int ex = 0; ex < Nex; ++ex){
    for(int site = 0; site < Nvol; ++site){

      iQ1 = matrix(iQ,site,ex);
  
      iQ2 = iQ1 * iQ1;
      iQ3 = iQ1 * iQ2;
   
      c0 = 0.0;
      c1 = 0.0;
      for(int cc = 0; cc < NC_; ++cc){
	c0 += iQ3.i(cc,cc);
	c1 += iQ2.r(cc,cc);
      }
      c0 = -c0/3.0;
      c1 = -c1/2.0;
      c0max = 2.0 * pow(c1/3.0,1.5);
  

      theta = acos(c0/c0max);
      u_val = sqrt(c1/3.0) * cos(theta/3.0);
      w = sqrt(c1) * sin(theta/3.0);
      xi0 = func_xi0(w);
      u2 = u_val * u_val;
      w2 = w * w;
      cosw = cos(w);

      ixi0 = dcomplex(0.0,xi0);
      emiu = dcomplex(cos(u_val),-sin(u_val));
      e2iu = dcomplex(cos(2.0*u_val),sin(2.0*u_val));

      h0 = e2iu * dcomplex(u2-w2,0.0)
	+ emiu * (dcomplex(8.0*u2*cosw,0.0)
		  + dcomplex(2.0*u_val*(3.0*u2 + w2),0.0)*ixi0);
      h1 = dcomplex(2*u_val,0.0) * e2iu
	- emiu*(dcomplex(2.0*u_val*cosw,0.0)
		- dcomplex(3.0*u2-w2,0.0) * ixi0);
      h2 = e2iu - emiu * (dcomplex(cosw,0.0)
			  + dcomplex(3.0*u_val,0.0)*ixi0);

      fden = 1.0/(9.0*u2 - w2);
      f0 = h0 * dcomplex(fden,0.0);
      f1 = h1 * dcomplex(fden,0.0);
      f2 = h2 * dcomplex(fden,0.0);

      for(int cc = 0; cc < NC_*NC_; ++cc){
	qt =  f0 * dcomplex(iQ0.r(cc), iQ0.i(cc))
	  + f1 * dcomplex(iQ1.i(cc),-iQ1.r(cc))
	  - f2 * dcomplex(iQ2.r(cc), iQ2.i(cc));
	e_iQ.data.set(e_iQ.format.index(2*cc,site,ex),  qt.real());
	e_iQ.data.set(e_iQ.format.index(2*cc+1,site,ex),qt.imag());
      }
      
    }
  }
}

//====================================================================
double Smear_Stout::func_xi0(double w) const{

  double xi0 = sin(w)/w;
  if( w < 1e-4 )
    CCIO::cout << "[Smear_stout] w too small: "<< w <<"\n";

  return xi0;

}

//====================================================================
//============================================================END=====
