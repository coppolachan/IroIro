/*! 
  @file dirac_clover.cpp
  @brief Dirac Clover operator class. Functions definition
*/
#include "dirac_clover.hpp"
#include "Communicator/comm_io.hpp"
#include "Tools/sunMatUtils.hpp"
#include "include/common_fields.hpp"

using namespace std;
using namespace FieldUtils;

typedef ShiftField_up<GaugeFieldFormat> FieldUP;//shifter up   for matrices
typedef ShiftField_dn<GaugeFieldFormat> FieldDN;//shifter down for matrices

void (Dirac_Clover::*Dirac_Clover::isigma[])(Field&,const Field&)const
= {&Dirac_Clover::isigma_diag,
   &Dirac_Clover::isigma_12,
   &Dirac_Clover::isigma_13,
   &Dirac_Clover::isigma_14,
   &Dirac_Clover::isigma_21,
   &Dirac_Clover::isigma_diag,
   &Dirac_Clover::isigma_23,
   &Dirac_Clover::isigma_24,
   &Dirac_Clover::isigma_31,
   &Dirac_Clover::isigma_32,
   &Dirac_Clover::isigma_diag,
   &Dirac_Clover::isigma_34,
   &Dirac_Clover::isigma_41,
   &Dirac_Clover::isigma_42,
   &Dirac_Clover::isigma_43,
   &Dirac_Clover::isigma_diag,};

void Dirac_Clover::isigma_diag(Field&,const Field&) const{
    CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::isigma.\n";
    abort();
}

void Dirac_Clover::mult_sw(Field& v_out, const Field& w) const {
  using namespace SUNmat_utils;
  using namespace SUNvec_utils;

  //Temporary variable declarations
  SUNvec v1;
  Field wt(fsize_);

  //Actual calculation - START
  isigma_23(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(Bx_,site) * v(wt,s,site);
      v_out.set(ff_.cslice(s,site), v1.getva());
    }
  }
  isigma_31(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(By_,site) * v(wt,s,site);
      v_out.add(ff_.cslice(s,site), v1.getva());
    }
  }
  isigma_12(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(Bz_,site) * v(wt,s,site);
      v_out.add(ff_.cslice(s,site),v1.getva());
    }
  }
  isigma_41(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(Ex_,site) * v(wt,s,site);
      v_out.add(ff_.cslice(s,site), v1.getva());
    }
  }
  isigma_42(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(Ey_,site) * v(wt,s,site);
      v_out.add(ff_.cslice(s,site),v1.getva());
    }
  }
  isigma_43(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(Ez_,site) * v(wt,s,site);
      v_out.add(ff_.cslice(s,site), v1.getva());
    }
  }
  v_out *= Dw->getKappa() * csw_;
}

//====================================================================
void Dirac_Clover::set_csw() {
  _Message(DEBUG_VERB_LEVEL, "[DiracClover] Setting Field Strenght\n");
  set_fieldstrength(Bx_,1,2);
  set_fieldstrength(By_,2,0);
  set_fieldstrength(Bz_,0,1);
  set_fieldstrength(Ex_,3,0);
  set_fieldstrength(Ey_,3,1);
  set_fieldstrength(Ez_,3,2);
}

//====================================================================
/*! @brief Calculates the term \f$F_{\mu.\nu}\f$  */
void Dirac_Clover::set_fieldstrength(GaugeField1D& field_strength,
				     const int mu, const int nu){
  using namespace SUNmat_utils;

  //.................. Temporary variables declaration
  GaugeField1D Cup, Cdn;
  GaugeField1D U_mu((*u_)[gf_.dir_slice(mu)]);  /*< @brief \f$U_\mu(x)\f$ */
  GaugeField1D w1, w2, v1, v2;
  //.......................................

  Cup.U = stpl_.upper(*u_,mu,nu); // Upper staple V_+mu
  Cdn.U = stpl_.lower(*u_,mu,nu); // Lower staple V_-mu

  for(int site = 0; site < Nvol_; ++site){
    w1.set_matrix(site,u(U_mu,site)   *u_dag(Cup,site));// U_mu(x)*(V_+mu)^dag
    w2.set_matrix(site,u(U_mu,site)   *u_dag(Cdn,site));// U_mu(x)*(V_-mu)^dag
    v1.set_matrix(site,u_dag(Cup,site)*u(U_mu,site)   );// (V_+mu)^dag*U_mu(x)
    v2.set_matrix(site,u_dag(Cdn,site)*u(U_mu,site)   );// (V_-mu)^dag*U_mu(x)
  }
  w1.U -= w2.U;
  //    +--<--+ 
  //    |     | w1
  // (x)+-->--+
  //    |     | w2
  //    +--<--+  

  v1.U -= v2.U; 
 
  FieldDN DnMu(v1, mu);
  //    +--<--+ 
  //    |     | v1
  //    +-->--+(x)
  //    |     | v2
  //    +--<--+  

  //Sum up the four terms
  w1.U += DnMu.getva();

  for(int site = 0; site < Nvol_; ++site)
    field_strength.U.set(field_strength.Format.cslice(0,site,0),
			 anti_hermite(u(w1,site)));
  field_strength.U *= 0.25;
}

/////////////////////////////////////////////////////////////////////

/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 1, \nu = 2\f$ */
void Dirac_Clover::isigma_12(Field& w,const Field& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site),-v[ff_.index_i(c,0,site)]);
      w.set(ff_.index_i(c,0,site), v[ff_.index_r(c,0,site)]);
      w.set(ff_.index_r(c,1,site), v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_i(c,1,site),-v[ff_.index_r(c,1,site)]);
      w.set(ff_.index_r(c,2,site),-v[ff_.index_i(c,2,site)]);
      w.set(ff_.index_i(c,2,site), v[ff_.index_r(c,2,site)]);
      w.set(ff_.index_r(c,3,site), v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_i(c,3,site),-v[ff_.index_r(c,3,site)]);
    }
  }
}
void Dirac_Clover::isigma_13(Field& w,const Field& v) const{
  isigma_31(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_14(Field& w,const Field& v) const{
  isigma_41(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_21(Field& w,const Field& v) const{
  isigma_12(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_23(Field& w,const Field& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site),-v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_i(c,0,site), v[ff_.index_r(c,1,site)]);
      w.set(ff_.index_r(c,1,site),-v[ff_.index_i(c,0,site)]);
      w.set(ff_.index_i(c,1,site), v[ff_.index_r(c,0,site)]);
      w.set(ff_.index_r(c,2,site),-v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_i(c,2,site), v[ff_.index_r(c,3,site)]);
      w.set(ff_.index_r(c,3,site),-v[ff_.index_i(c,2,site)]);
      w.set(ff_.index_i(c,3,site), v[ff_.index_r(c,2,site)]);
    }
  }
}
void Dirac_Clover::isigma_24(Field& w,const Field& v) const{
  isigma_42(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_31(Field& w,const Field& v) const{

  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site), v[ff_.index_r(c,1,site)]);
      w.set(ff_.index_i(c,0,site), v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_r(c,1,site),-v[ff_.index_r(c,0,site)]);
      w.set(ff_.index_i(c,1,site),-v[ff_.index_i(c,0,site)]);
      w.set(ff_.index_r(c,2,site), v[ff_.index_r(c,3,site)]);
      w.set(ff_.index_i(c,2,site), v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_r(c,3,site),-v[ff_.index_r(c,2,site)]);
      w.set(ff_.index_i(c,3,site),-v[ff_.index_i(c,2,site)]);
    }
  }
}
void Dirac_Clover::isigma_32(Field& w,const Field& v) const{
  isigma_23(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_34(Field& w,const Field& v) const{
  isigma_43(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_41(Field& w,const Field& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site), v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_i(c,0,site),-v[ff_.index_r(c,3,site)]);
      w.set(ff_.index_r(c,1,site), v[ff_.index_i(c,2,site)]);
      w.set(ff_.index_i(c,1,site),-v[ff_.index_r(c,2,site)]);
      w.set(ff_.index_r(c,2,site), v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_i(c,2,site),-v[ff_.index_r(c,1,site)]);
      w.set(ff_.index_r(c,3,site), v[ff_.index_i(c,0,site)]);
      w.set(ff_.index_i(c,3,site),-v[ff_.index_r(c,0,site)]);
    }
  }
}
void Dirac_Clover::isigma_42(Field& w,const Field& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site),-v[ff_.index_r(c,3,site)]);
      w.set(ff_.index_i(c,0,site),-v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_r(c,1,site), v[ff_.index_r(c,2,site)]);
      w.set(ff_.index_i(c,1,site), v[ff_.index_i(c,2,site)]);
      w.set(ff_.index_r(c,2,site),-v[ff_.index_r(c,1,site)]);
      w.set(ff_.index_i(c,2,site),-v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_r(c,3,site), v[ff_.index_r(c,0,site)]);
      w.set(ff_.index_i(c,3,site), v[ff_.index_i(c,0,site)]);
    }
  }
}
void Dirac_Clover::isigma_43(Field& w,const Field& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.set(ff_.index_r(c,0,site), v[ff_.index_i(c,2,site)]);
      w.set(ff_.index_i(c,0,site),-v[ff_.index_r(c,2,site)]);
      w.set(ff_.index_r(c,1,site),-v[ff_.index_i(c,3,site)]);
      w.set(ff_.index_i(c,1,site), v[ff_.index_r(c,3,site)]);
      w.set(ff_.index_r(c,2,site), v[ff_.index_i(c,0,site)]);
      w.set(ff_.index_i(c,2,site),-v[ff_.index_r(c,0,site)]);
      w.set(ff_.index_r(c,3,site),-v[ff_.index_i(c,1,site)]);
      w.set(ff_.index_i(c,3,site), v[ff_.index_r(c,1,site)]);
    }
  }
}

const Field Dirac_Clover::gamma5(const Field& f) const{
  return Dw->gamma5(f);
}

const Field Dirac_Clover::mult(const Field& f) const{
  Field w (ff_.size());
  Field w2(ff_.size());

  w  = Dw->mult(f);
  mult_sw(w2,f);
  w -= w2;
  return w;
}

const Field Dirac_Clover::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

//====================================================================
const Field Dirac_Clover::md_force(const Field& eta,const Field& zeta)const{
  //Wilson term
  Field force = Dw->md_force(eta,zeta);
 
  force += md_force_block(eta, zeta);
  force += md_force_block(zeta, eta);

  return force;
}
//====================================================================
const Field Dirac_Clover::md_force_block(const Field& eta,
					 const Field& zeta)const{
  using namespace SUNmat_utils;
  using namespace SUNvec_utils;

  //.................... Temporary variables declaration
  int Nd = CommonPrms::instance()->Nd();

  GaugeField force;
  GaugeField1D fce_tmp1, fce_tmp2;
  Field eta2(fsize_), vleft(fsize_), vright(fsize_);
  Field shifted(fsize_);
  GaugeField1D U_mu, U_nu;
  GaugeField1D Cmu_up, Cnu_up;
  SUNvec vect;
  //...................................................
  
  //Clover term: 8 terms, see Matsufuru-san's note

  for(int mu = 0; mu < Ndim_; ++mu){
    for(int nu = 0; nu < Ndim_; ++nu){
      if(nu == mu) continue;

      (this->*isigma[Ndim_*mu+nu])(eta2,eta); //sigma_{mu,nu} mult

      U_nu.U = (*u_)[gf_.dir_slice(nu)];
      U_mu.U = (*u_)[gf_.dir_slice(mu)]; 

      //term 1 and 5 --------------------------------------
      Cmu_up.U  = stpl_.upper(*u_, mu,nu); //V_mu
      Cmu_up.U -= stpl_.lower(*u_, mu,nu); //V_+mu - V_-mu 
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = u_dag(Cmu_up,site) * v(eta2,s,site);
	  vect = u(U_mu,site) * vect;
	  vright.set(ff_.cslice(s,site), vect.getva());
	} 
      }
      fce_tmp2.U = field_oprod(zeta,vright); 

      //term 2 --------------------------------------
      Cnu_up.U = stpl_.upper(*u_, nu, mu);
      sf_up_[nu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(Cnu_up.U,site))* v(sf_up_[nu],s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      }
      sf_up_[nu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      }// U_nu(x)*zeta(x+nu)
      fce_tmp2.U += field_oprod(vleft,vright);

      //term 4 and 8 --------------------------------------
      sf_up_[mu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      } // U_mu (x) * eta2(x+mu)   
      sf_up_[mu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(Cmu_up.U,site))* v(sf_up_[mu],s,site);
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      } //  (V_+mu + V_-mu) * zeta(x+mu)
      fce_tmp2.U += field_oprod(vleft,vright);
      
      //term 3 --------------------------------------
      sf_up_[nu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      }       
      sf_up_[mu]->setf(vright); //U_nu(x+mu)*eta2(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      } //vright = U_mu(x)*U_nu(x+mu)*eta2(x+mu+nu)  
     
      sf_up_[mu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      }        
      sf_up_[nu]->setf(vleft); //U_mu(x+nu)*zeta(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      } //vright = U_nu(x)*U_mu(x+nu)*zeta(x+mu+nu)  
      fce_tmp2.U += field_oprod(vleft,vright);

      //term 6 --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(Cnu_up.U,site))* v(eta2,s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      }       
      sf_dn_[nu]->setf(vright); // V_+nu (x-nu)*eta2(x-nu)
      shifted = sf_dn_[nu]->getva();

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(U_nu.U,site))* v(zeta,s,site);
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      }   
      sf_dn_[nu]->setf(vleft);
      fce_tmp2.U -= field_oprod((Field)(sf_dn_[nu]->getva()), shifted);

      // term 7  --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(U_nu.U,site))* v(eta2,s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      }   
      sf_dn_[nu]->setf(vright);
      shifted = sf_dn_[nu]->getva(); // shifted = Udag_nu(x-nu)*eta2(x-nu)
      sf_up_[mu]->setf(shifted);// Udag_nu(x+mu-nu)*eta2(x+mu-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_.cslice(s,site), vect.getva());
	}
      } //vright = U_mu(x)*Udag_nu(x+mu-nu)*eta2(x+mu-nu)
     
      sf_up_[mu]->setf(zeta); //zeta(x+mu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vect = (u_dag(U_nu.U,site))*vect;
	  vleft.set(ff_.cslice(s,site), vect.getva());
	}
      }       
      sf_dn_[nu]->setf(vleft);//sf_dn_[nu] = Udag_nu(x-nu)*U_mu(x-nu)*zeta(x+mu-nu)
      fce_tmp2.U -= field_oprod((Field)(sf_dn_[nu]->getva()), vright);
      fce_tmp2.U *= - Dw->getKappa() * csw_ / 8.0; 

      SUNmat force_mat;
      for(int site = 0; site<Nvol_; ++site){
	force_mat = u(fce_tmp2.U,site);
	force.U.add(force.Format.cslice(0,site,mu), force_mat.getva() ); 
      }
    }
  }
  return force.U;
}
//====================================================================
const vector<int> Dirac_Clover::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}
