/*!
  @file action_gauge_rect.cpp

  @brief Definition of the ActionGaugeRect class
*/
#include "action_gauge_rect.hpp"
#include "Communicator/comm_io.hpp"

using namespace std;

typedef ShiftField_up<GaugeFieldFormat> FieldUP;
typedef ShiftField_dn<GaugeFieldFormat> FieldDN;

double ActionGaugeRect::calc_H(){
  using namespace SUNmat_utils;
  int NP = CommonPrms::instance()->NP();
  double plaqF = 0.0;
  double rectF = 0.0;

  GaugeField1D Cup1, Cup2;
  valarray<double> w(sf_->size()), v(sf_->size()), res(sf_->size());

  // From Matsufuru-san code
  for(int mu = 0; mu < Ndim_; ++mu){
    FieldUP UpMu(sf_, mu);//shifter Up in mu direction

   for(int nu = mu+1; nu < Ndim_; ++nu){  
     FieldUP UpNu(sf_, nu);//shifter Up in nu direction
     Cup1.U = stpl_->upper(*u_,mu,nu);
     Cup2.U = stpl_->upper(*u_,nu,mu);
     
     // plaquette term
     for(int site = 0; site < Nvol_; ++site){
        plaqF += ReTr( u(*u_,site,mu) * u_dag(Cup1,site) );
     }

     // rectangular terms
     //       mu,v (UpNu)                               
     //      +-->--+-->--+                                                    
     // nu,w |           |Cup2_dag(site+mu) (UpMu)
     //  site+     +--<--+                               
 
     w = (*u_)[gf_.dir_slice(nu)];
     v = (*u_)[gf_.dir_slice(mu)]; 
     UpNu.setf(v);  
     UpMu.setf(Cup2.U);      

     res = 0.0;
     for(int site = 0; site<Nvol_; ++site){
       res[sf_->cslice(0,site)] = (u(w,*sf_,site)*u(UpNu,site)*u_dag(UpMu,site)).getva();
       rectF += ReTr(u(*u_,gf_,site,mu)* u_dag(res,site));
     }

     //     Cup1(site+nu)
     //       +-->--+
     //       |     |
     //       +     +
     // nu,w  |     | w+mu (UpMu)
     // site  +     +

     w = (*u_)[gf_.dir_slice(nu)];
     UpNu.setf(Cup1.U);
     UpMu.setf(w);
  
     res = 0.0;
     for(int site = 0; site<Nvol_; ++site){
       res[sf_->cslice(0,site)] = (u(w,*sf_,site)*u(UpNu,site)*u_dag(UpMu,site)).getva();
       rectF += ReTr(u(*u_,gf_,site,mu)* u_dag(res,site));
     }     
   }
  }


  plaqF = Communicator::instance()->reduce_sum(plaqF);
  rectF = Communicator::instance()->reduce_sum(rectF);

  double plaq = plaqF/Nc_;
  double rect = rectF/Nc_;
  
  CCIO::cout<<" -- Plaquette = "<< plaqF/(Nvol_*NP) << "\n";
  
  double Hgauge = Params.c_plaq*(Nvol_*NP*Ndim_*(Ndim_-1.0)/2.0 - plaq)
    + Params.c_rect*(Nvol_*NP*Ndim_*(Ndim_-1.0) - rect);
    
    
#if VERBOSITY>=ACTION_VERB_LEVEL
    CCIO::cout << "[ActionGaugeRect] H = "<< Hgauge << "\n";
#endif
  
  return Hgauge;
}

Field ActionGaugeRect::md_force(const void*){
  using namespace SUNmat_utils;
  
  Field force(gf_.size());
  GaugeField1D force_pl;
  GaugeField1D force_rect;

  GaugeField1D Cup1;
  GaugeField1D Cup2;
  GaugeField1D Cdn1;
  GaugeField1D Cdn2;
  //valarray for speed purposes (GaugeField is slower for these)
  valarray<double> w(sf_->size()), v(sf_->size()), res(sf_->size());
  for(int mu = 0; mu < Ndim_; ++mu){
    force_pl.U   = 0.0;
    force_rect.U = 0.0;
    FieldUP UpMu(sf_, mu);//shifter Up in mu direction
    for(int nu = 0; nu < Ndim_; ++nu){
      if (nu == mu) continue;
      FieldUP UpNu(sf_, nu);//shifter Up in nu direction
      FieldDN DnNu(sf_, nu);//shifter Down in nu direction

      Cup1.U = stpl_->upper(*u_,mu,nu);
      Cup2.U = stpl_->upper(*u_,nu,mu);
      Cdn1.U = stpl_->lower(*u_,mu,nu);
      Cdn2.U = stpl_->lower(*u_,nu,mu);

      // plaquette term
      force_pl.U += Cup1.U;
      force_pl.U += Cup2.U;
      force_pl.U *= Params.c_plaq;
      // rectangular terms
      // ^nu
      // |  
      // +-->mu
      //
      // (x) is the site position
      
      w = (*u_)[gf_.dir_slice(mu)];  //U_mu links
      v = (*u_)[gf_.dir_slice(nu)];  //U_nu links    

      //        w
      //      +-->--+-->--+
      //   v  |           |   term  (Cup2)
      //     (x)    +--<--+      
  
      UpMu.setf(Cup2.U); // Cup2(x+mu)
      UpNu.setf(w);      // U_mu(x+nu)

      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u(v,*sf_,site)*u(UpNu,site)*u_dag(UpMu,site)).getva();
      }       
      force_rect.U += res;

      //      +-->--+
      //      |     |
      //      +     +   term
      //   v  |     |  v(x+mu)
      //     (x)    v

      UpMu.setf(v);    // v(x+mu)
      UpNu.setf(Cup1.U); // Cup1(x+nu)
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u(v,*sf_,site)*u(UpNu,site)*u_dag(UpMu,site)).getva();
      }       
      force_rect.U += res;     
    
      //      +-->--+-->--+
      //      |           |   term
      //      +--<-(x)    v
      
      UpNu.setf(w); // w(x+nu)
  
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u(Cdn2.U,site)*u(UpNu,site)*u_dag(UpMu,site)).getva();
      }       
      force_rect.U += res;

      //     (x)    +--<--+
      //      |           |   term
      //      +-->--+-->--+

      UpMu.setf(Cup2.U);
      
   
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u_dag(v,*sf_,site)*u(w,*sf_,site)*u(UpMu,site)).getva();
      } 
      DnNu.setf(res); //res(x-nu)
      force_rect.U += DnNu.getva();

      //     (x)    ^
      //      |     |
      //      +     +   term
      //      |     |
      //      +-->--+

      UpMu.setf(v); //v(x+mu)

      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u_dag(v,*sf_,site)*u(Cdn1,site)*u(UpMu,site)).getva();
      } 
      DnNu.setf(res); //res(x-nu)     
      force_rect.U += DnNu.getva();
  
      //      +--<-(x)    ^
      //      |           |   term
      //      +--<--+--<--+
      
      res = 0.0;
      for(int site = 0; site<Nvol_; ++site){
	res[sf_->cslice(0,site)] = (u_dag(Cdn2,site)*u(w,*sf_,site)*u(UpMu,site)).getva();
      }       
      DnNu.setf(res); //res(x-nu)     
      force_rect.U += DnNu.getva();

      force_rect.U *= Params.c_rect;
      force_rect.U += force_pl.U; //total force (staples term)
    }

    
    for(int site = 0; site<Nvol_; ++site){
      force[sf_->cslice(0,site)] = (u(*u_,gf_,site)*u_dag(force_rect.U,site)).getva();
    }   


  force *= -Params.beta/Nc_;
  }

#if VERBOSITY>=ACTION_VERB_LEVEL
  monitor_force(force, "ActionGaugeRect");
#endif

  return force;
}
