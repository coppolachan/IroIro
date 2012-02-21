/*!
  @file action_gauge_rect.cpp

  @brief Definition of the ActionGaugeRect class
*/
#include "action_gauge_rect.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "lib/Main/Geometry/mapper.hpp"
#include "include/messages_macros.hpp"

using namespace SUNmatUtils;
using namespace FieldUtils;


double ActionGaugeRect::calc_H(){
  using namespace MapsEnv;

  int NP = CommonPrms::instance()->NP();
  double plaqF = 0.0;
  double rectF = 0.0;

  const Staples stpl_;
  GaugeField1D Cup1, Cup2;
  GaugeField1D U_nu, U_mu, res;
  GaugeField1D UpNu, UpMu;

  // From Matsufuru-san code
  for(int mu = 0; mu < NDIM_; ++mu){

   for(int nu = mu+1; nu < NDIM_; ++nu){  
     Cup1 = stpl_.upper(*u_,mu,nu);
     Cup2 = stpl_.upper(*u_,nu,mu);
     
     // plaquette term
     for(int site = 0; site < Nvol_; ++site){
       plaqF += ReTr( matrix(*u_,site,mu) * matrix_dag(Cup1,site) );
     }
     
     // rectangular terms
     //       mu,U_mu (UpNu) 
     //         +-->--+-->--+
     // nu,U_nu |           |Cup2_dag(site+mu) (UpMu)
     //     site+     +--<--+                               
 
     U_nu = DirSlice(*u_, nu);
     U_mu = DirSlice(*u_, mu); 
     UpNu = shift(U_mu, nu, Forward);      //U_mu(x+nu)
     UpMu = shift(Cup2, mu, Forward);      //Cup2(x+mu)     

     res = 0.0;
     for(int site = 0; site<Nvol_; ++site){
       res.data[res.format.cslice(0,site)] = 
	 (matrix(U_nu,site) * matrix(UpNu,site) * matrix_dag(UpMu,site)).getva();
       rectF += ReTr(matrix(*u_,site,mu)* matrix_dag(res,site));
     }

     //        Cup1(site+nu)
     //          +-->--+
     //          |     |
     //          +     +
     // nu,U_nu  |     | U_nu+mu (UpMu)
     //    site  +     +

     UpNu = shift(Cup1, nu, Forward);   //Cup1(x+nu)
     UpMu = shift(U_nu, mu, Forward);   //U_nu(x+mu)
  
     res = 0.0;
     for(int site = 0; site<Nvol_; ++site){
       res.data[res.format.cslice(0,site)] =
	 (matrix(U_nu,site) * matrix(UpNu,site) * matrix_dag(UpMu,site)).getva();
       rectF += ReTr(matrix(*u_,site,mu)* matrix_dag(res,site));
     } 
         
   }
  }

  plaqF = Communicator::instance()->reduce_sum(plaqF);
  rectF = Communicator::instance()->reduce_sum(rectF);

  double plaq = plaqF/NC_;
  double rect = rectF/NC_;

  double Hgauge = Params.c_plaq*(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0 - plaq)
                + Params.c_rect*(Nvol_*NP*NDIM_*(NDIM_-1.0)     - rect);

  Hgauge *= Params.beta;
    
  _Message(ACTION_VERB_LEVEL, "    [ActionGaugeRect] H = "<<Hgauge<<"\n");
  _Message(1,"    -- Plaquette = "<< plaq/(Nvol_*NP*NDIM_*(NDIM_-1.0)/2.0) << "\n");
  return Hgauge;
}

GaugeField ActionGaugeRect::md_force(){
 using namespace MapsEnv;
  
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
  GaugeField1D U_mu, U_nu, res;
  GaugeField1D UpMu, UpNu;

  for(int mu = 0; mu < NDIM_; ++mu){
    force_pl   = 0.0;
    force_rect = 0.0;
    //FieldUP UpMu(sf_, mu);//shifter Up in mu direction
    for(int nu = 0; nu < NDIM_; ++nu){
      if (nu == mu) continue;
      //FieldUP UpNu(sf_, nu);//shifter Up in nu direction
      //FieldDN DnNu(sf_, nu);//shifter Down in nu direction

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
      
      U_mu = DirSlice(*u_, mu);   //U_mu links
      U_nu = DirSlice(*u_, nu);   //U_nu links    

      //          U_mu
      //         +-->--+-->--+
      //   U_nu  |           |   term  (Cup2)
      //        (x)    +--<--+      
  
      UpMu = shift(Cup2, mu, Forward); // Cup2(x+mu)
      UpNu = shift(U_mu, nu, Forward); // U_mu(x+nu)

      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (matrix(U_nu,site)*matrix(UpNu,site)*matrix_dag(UpMu,site)).getva();
      }       
      force_rect += res;

      //         +-->--+
      //         |     |
      //         +     +   term
      //   U_nu  |     |  U_nu(x+mu) (UpMu)
      //        (x)    v

      UpMu = shift(U_nu, mu, Forward);    // U_nu(x+mu)
      UpNu = shift(Cup1, nu, Forward);    // Cup1(x+nu)
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] =
	  (matrix(U_nu,site)*matrix(UpNu,site)*matrix_dag(UpMu,site)).getva();
      }       
      force_rect += res;     
    
      //           U_mu(x+nu)
      //      +-->--+-->--+
      //      |           |   term
      //      +--<-(x)    v
      
      UpNu = shift(U_mu, nu, Forward); // U_mu(x+nu)
  
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (matrix(Cdn2,site)*matrix(UpNu,site)*matrix_dag(UpMu,site)).getva();
      }       
      force_rect += res;

      //     (x)    +--<--+
      //      |           |   term
      //      +-->--+-->--+

      UpMu = shift(Cup2, mu, Forward);

      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (matrix_dag(U_nu,site)*matrix(U_mu,site)*matrix(UpMu,site)).getva();
      } 
      force_rect += shift(res, nu, Backward); //+=res(x-nu)
    

      //     (x)    ^
      //      |     |
      //      +     +   term
      //      |     |
      //      +-->--+

      UpMu = shift(U_nu, mu, Forward); //U_nu(x+mu)

      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (matrix_dag(U_nu,site)*matrix(Cdn1,site)*matrix(UpMu,site)).getva();
      } 
      force_rect += shift(res, nu, Backward);//+=res(x-nu)    
  
      //      +--<-(x)    ^
      //      |           |   term
      //      +--<--+--<--+
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res.data[res.format.cslice(0,site)] = 
	  (matrix_dag(Cdn2,site)*matrix(U_mu,site)*matrix(UpMu,site)).getva();
      }       
      force_rect += shift(res, nu, Backward);//+=res(x-nu)  

    }

    force_pl   *= Params.c_plaq;
    force_rect *= Params.c_rect;
    force_rect += force_pl; //force_rect = total force (staples term)
    for(int site = 0; site<Nvol_; ++site){
      force_mat = (matrix(*u_,site,mu)*matrix_dag(force_rect,site));
      SetMatrix(force, anti_hermite(force_mat), site, mu);
    } 
  }

  force *= 0.5*Params.beta/NC_;

  _MonitorMsg(ACTION_VERB_LEVEL, Action, force, "ActionGaugeRect");

  return force;
}
