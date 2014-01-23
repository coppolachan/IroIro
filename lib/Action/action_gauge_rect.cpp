/*!
  @file action_gauge_rect.cpp
  @brief Definition of the ActionGaugeRect class

  Time-stamp: <2014-01-15 14:09:29 cossu>
*/

#include "action_gauge_rect.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/messages_macros.hpp"
#include "Measurements/GaugeM/polyakovLoop.hpp"

using namespace SUNmatUtils;
using namespace FieldUtils;



double ActionGaugeRect::calc_H(){
  using namespace Mapping;

  const Staples stpl_;
  GaugeField1D Cup1, Cup2;
  GaugeField1D res;
  GaugeField1D UpNu, UpMu;
  GaugeField1D U_mu, U_nu;

  double plaqF = 0.0;
  double rectF = 0.0;

  for(int mu=0; mu<NDIM_; ++mu){
    U_mu = DirSlice(*u_,mu); 
    
    for(int nu = mu+1; nu < NDIM_; ++nu){
      Cup1 = stpl_.upper(*u_,mu,nu);
      Cup2 = stpl_.upper(*u_,nu,mu);
      
      // plaquette term
      for(int site=0; site<Nvol_; ++site)
	plaqF += ReTr(mat(*u_,site,mu)*mat_dag(Cup1,site));
      
      U_nu = DirSlice(*u_,nu);
      // rectangular terms
      //       mu,U_mu (UpNu) 
      //         +-->--+-->--+
      // nu,U_nu |           |Cup2_dag(site+mu) (UpMu)
      //     site+     +--<--+                               
      
      UpNu = shiftField(U_mu,nu,Forward());      //U_mu(x+nu)
      UpMu = shiftField(Cup2,mu,Forward());      //Cup2(x+mu)     
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
	rectF += ReTr(mat(*u_,site,mu)*mat_dag(res,site));
      }
      
      //        Cup1(site+nu)
      //          +-->--+
      //          |     |
      //          +     +
      // nu,U_nu  |     | U_nu+mu (UpMu)
      //    site  +     +
      
      UpNu = shiftField(Cup1,nu,Forward());   //Cup1(x+nu)
      UpMu = shiftField(U_nu,mu,Forward());   //U_nu(x+mu)
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] =
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
	rectF += ReTr(mat(*u_,site,mu)*mat_dag(res,site));
      }
    }
  }

  Communicator::instance()->sync();
  plaqF = Communicator::instance()->reduce_sum(plaqF);
  rectF = Communicator::instance()->reduce_sum(rectF);
 
  int NP = CommonPrms::instance()->NP(); 
  double plaq = plaqF/NC_;
  double rect = rectF/NC_;
  double Hgauge = Params.c_plaq*(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0 -plaq)
    + Params.c_rect*(Nvol_*NP*NDIM_*(NDIM_-1.0) -rect);
  Hgauge *= Params.beta;
  
  _Message(ACTION_VERB_LEVEL, "    [ActionGaugeRect] H = "<<Hgauge<<"\n");
  _Message(BASE_VERB_LEVEL,   "    [ActionGaugeRect] Plaquette = "
	   << plaq/(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0) << "\n");

  //Measure other observables
  
  PolyakovLoop PLmeas(TDIR);
  std::complex<double> pl = PLmeas.calc_SUN(*u_);
  _Message(ACTION_VERB_LEVEL,"    [ActionGaugeRect] PL = "<<pl.real()<< "  "<< pl.imag() << "\n");

  return Hgauge;
}

GaugeField ActionGaugeRect::md_force(){
 using namespace Mapping;
  
  const Staples stpl_;
  SUNmat force_mat;
  GaugeField force;
  GaugeField1D force_pl;
  GaugeField1D force_rect;

  GaugeField1D Cup1;
  GaugeField1D Cup2;
  GaugeField1D Cdn1;
  GaugeField1D Cdn2;

  //check speed
  GaugeField1D res;
  GaugeField1D UpMu, UpNu;
  GaugeField1D U_mu, U_nu;

  for(int mu=0; mu<NDIM_; ++mu){
    force_pl   = 0.0;
    force_rect = 0.0;

    U_mu = DirSlice(*u_,mu);   //U_mu links

    for(int nu=0; nu<NDIM_; ++nu){
      if (nu == mu) continue;
      
      Cup1 = stpl_.upper(*u_,mu,nu);
      Cup2 = stpl_.upper(*u_,nu,mu);
      Cdn1 = stpl_.lower(*u_,mu,nu);
      Cdn2 = stpl_.lower(*u_,nu,mu);

      // plaquette term
      force_pl += Cup1;
      force_pl += Cdn1;

      // rectangular terms
      // ^nu
      // |  
      // +-->mu
      //
      // (x) is the site position
      U_nu = DirSlice(*u_,nu);   //U_nu links    
      //          U_mu
      //         +-->--+-->--+
      //   U_nu  |           |   term  (Cup2)
      //        (x)    +--<--+      

      UpMu = shiftField(Cup2,mu,Forward()); // Cup2(x+mu)
      UpNu = shiftField(U_mu,nu,Forward()); // U_mu(x+nu)

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
        res.data[res.format.cslice(0,site)] = 
          (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();

      force_rect += res;  
  
      //         +-->--+
      //         |     |
      //         +     +   term
      //   U_nu  |     |  U_nu(x+mu) (UpMu)
      //        (x)    v
      UpMu = shiftField(U_nu, mu, Forward());    // U_nu(x+mu)
      UpNu = shiftField(Cup1, nu, Forward());    // Cup1(x+nu)
      res = 0.0;

      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] =
	  (mat(U_nu,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();
      
     force_rect += res;     
      //           U_mu(x+nu)
      //      +-->--+-->--+
      //      |           |   term
      //      +--<-(x)    v
      UpNu = shiftField(U_mu,nu,Forward()); // U_mu(x+nu)

      res = 0.0; 
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat(Cdn2,site)*mat(UpNu,site)*mat_dag(UpMu,site)).getva();

      force_rect += res;
      //     (x)    +--<--+
      //      |           |   term
      //      +-->--+-->--+
      UpMu = shiftField(Cup2,mu,Forward());

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(U_nu,site)*mat(U_mu,site)*mat(UpMu,site)).getva();

      force_rect += shiftField(res,nu,Backward()); //+=res(x-nu)
      //     (x)    ^
      //      |     |
      //      +     +   term
      //      |     |
      //      +-->--+
      UpMu = shiftField(U_nu, mu, Forward()); //U_nu(x+mu)

      res = 0.0;
      for(int site=0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(U_nu,site)*mat(Cdn1,site)*mat(UpMu,site)).getva();
      } 
      force_rect += shiftField(res,nu,Backward());//+=res(x-nu)    
      //      +--<-(x)    ^
      //      |           |   term
      //      +--<--+--<--+

      res = 0.0;
      for(int site=0; site<Nvol_; ++site)
	res.data[res.format.cslice(0,site)] = 
	  (mat_dag(Cdn2,site)*mat(U_mu,site)*mat(UpMu,site)).getva();

      force_rect += shiftField(res,nu,Backward());//+=res(x-nu)  
    }
    force_pl   *= Params.c_plaq;
    force_rect *= Params.c_rect;
    force_rect += force_pl; //force_rect = total force (staples term)

    for(int site=0; site<Nvol_; ++site){
      force_mat = (mat(*u_,site,mu)*mat_dag(force_rect,site));
      SetMat(force, anti_hermite_traceless(force_mat), site, mu);
    }
  } 

  force *= 0.5*Params.beta/NC_;
  _MonitorMsg(ACTION_VERB_LEVEL,Action,force,"ActionGaugeRect");

  return force;
}


