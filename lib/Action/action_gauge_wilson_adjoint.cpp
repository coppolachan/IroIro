/*!
  @file action_gauge_wilson_adjoint.testing.cpp
  @brief Definition of the ActionGaugeWilsonAdjoint class (testing state)
*/
#include "action_gauge_wilson_adjoint.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/messages_macros.hpp"

double ActionGaugeWilsonAdjoint::calc_H(){
  //Number of plaquettes
  int Nplaq = CommonPrms::instance()->Lvol()*NDIM_*(NDIM_-1)/2.0;
  double plaq = stpl_.plaquette_adj(*u_);
  double Hgauge = Params.beta*Nplaq*(1.0-plaq);

  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilsonAdjoint] H = "<<Hgauge <<"\n");
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeWilsonAdjoint] Plaquette Adjoint = "<<plaq<<"\n");
  
  return Hgauge;
}

GaugeField ActionGaugeWilsonAdjoint::md_force(){
  using namespace FieldUtils;
  using namespace SUNmatUtils;
  using namespace Mapping;

  SUNmat pl, pl_sum, unitNc;
  GaugeField force;
  GaugeField1D tmp_up, tmp_dn; 
 
  GaugeField1D v, w, c;
  GaugeField1D WupMu, VupNu;

  unitNc = unity();
  unitNc /= (double)NC_;
  std::complex<double> trace;

  for(int m = 0; m < NDIM_; ++m){
    tmp_up = 0.0;
    tmp_dn = 0.0;
    for(int n=0; n< NDIM_; ++n) {
      if(n != m)
	{
	  tmp_up = stpl_.upper(*u_,m,n);
	  tmp_dn = stpl_.lower(*u_,m,n);
	  
	  for(int site=0; site < Nvol_; ++site){
	    pl_sum.zero();
	    pl = mat(*u_,site,m)*mat_dag(tmp_up,site);
	    //trace.real(ReTr(pl));
	    //trace.imag(-ImTr(pl));
	    trace = (ReTr(pl),-ImTr(pl));
	    pl_sum += pl*trace;
	    pl = mat(*u_,site,m)*mat_dag(tmp_dn,site);
	    //trace.real(ReTr(pl));
	    //trace.imag(ImTr(pl));
	    trace = (ReTr(pl),ImTr(pl));
	    pl_sum += pl*trace;

	    AddMat(force, anti_hermite_traceless(pl_sum), site, m);
	  }	
	}
    }
  }


  force *= Params.beta/(NC_*NC_-1.0);

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeWilsonAdjoint");

  return force;
}


