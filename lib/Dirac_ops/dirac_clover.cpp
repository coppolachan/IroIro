/*! 
  @file dirac_clover.cpp
  @brief Dirac Clover operator class. Functions definition
*/
#include "dirac_clover.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/common_fields.hpp"
#include "include/messages_macros.hpp"

using namespace std;

void (Dirac_Clover::*Dirac_Clover::isigma[])(FermionField&,
					     const FermionField&)const
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

void Dirac_Clover::isigma_diag(FermionField&,const FermionField&) const{
    CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::isigma.\n";
    abort();
}

void Dirac_Clover::mult_sw(FermionField& v_out, const FermionField& w) const {
  using namespace SUNmatUtils;
  using namespace SUNvec_utils;
  using namespace FieldUtils;

  //Temporary variable declarations
  SUNvec v1;
  FermionField wt;

  //Actual calculation - START
  isigma_23(wt,w);
  for(int site=0; site<Nvol_; ++site)
    for(int s=0; s<NDIM_; ++s)
      SetVec(v_out,(mat(Bx_,site)*vec(wt,s,site)),s,site);

  isigma_31(wt,w);
  for(int site=0; site<Nvol_; ++site)
    for(int s=0; s<NDIM_; ++s)
      AddVec(v_out,(mat(By_,site)*vec(wt,s,site)),s,site);
  
  isigma_12(wt,w);
  for(int site=0; site<Nvol_; ++site)
    for(int s=0; s<NDIM_; ++s)
      AddVec(v_out,(mat(Bz_,site)*vec(wt,s,site)),s,site);

  isigma_41(wt,w);
  for(int site=0; site<Nvol_; ++site)
    for(int s=0; s<NDIM_; ++s)
      AddVec(v_out,(mat(Ex_,site)*vec(wt,s,site)),s,site);

  isigma_42(wt,w);
  for(int site=0; site<Nvol_; ++site)
    for(int s=0; s<NDIM_; ++s)
      AddVec(v_out,(mat(Ey_,site)*vec(wt,s,site)),s,site);

  isigma_43(wt,w);
  for(int site=0; site<Nvol_; ++site){
    for(int s=0; s<NDIM_; ++s){
      v1 = mat(Ez_,site)*vec(wt,s,site);
      AddVec(v_out,v1,s,site);
    }
  }
  v_out *= Dw->getKappa() * csw_;
}

//====================================================================
void Dirac_Clover::set_csw() {
  _Message(DEBUG_VERB_LEVEL, "[DiracClover] Setting Field Strenght\n");

  set_fieldstrength(Bx_, 1, 2);
  set_fieldstrength(By_, 2, 0);
  set_fieldstrength(Bz_, 0, 1);
  set_fieldstrength(Ex_, 3, 0);
  set_fieldstrength(Ey_, 3, 1);
  set_fieldstrength(Ez_, 3, 2);
}
//====================================================================
/*! @brief Calculates the term \f$F_{\mu.\nu}\f$  */
void Dirac_Clover::set_fieldstrength(GaugeField1D& field_strength,
				     const int mu, const int nu){
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  using namespace Mapping;

  //.................. Temporary variables declaration
  GaugeField1D Cup, Cdn;
  GaugeField1D U_mu;  /*< @brief \f$U_\mu(x)\f$ */
  GaugeField1D w1, w2, v1, v2;
  Staples stpl;

  //really temporary - to be eliminated in future
  GaugeField temp_u_(*u_);
  //.......................................
  U_mu = DirSlice(temp_u_, mu);

  Cup = stpl.upper(temp_u_,mu,nu); // Upper staple V_+mu
  Cdn = stpl.lower(temp_u_,mu,nu); // Lower staple V_-mu

  for(int site = 0; site < Nvol_; ++site){
    SetMat(w1,mat(U_mu,site)   *mat_dag(Cup,site),site);// U_mu(x)*(V_+mu)^dag
    SetMat(w2,mat(U_mu,site)   *mat_dag(Cdn,site),site);// U_mu(x)*(V_-mu)^dag
    SetMat(v1,mat_dag(Cup,site)*mat(U_mu,site)   ,site);// (V_+mu)^dag*U_mu(x)
    SetMat(v2,mat_dag(Cdn,site)*mat(U_mu,site)   ,site);// (V_-mu)^dag*U_mu(x)
  }
  w1 -= w2;
  //    +--<--+ 
  //    |     | w1
  // (x)+-->--+
  //    |     | w2
  //    +--<--+  

  v1 -= v2; 
  //    +--<--+ 
  //    |     | v1
  //    +-->--+(x)
  //    |     | v2
  //    +--<--+  

  //Sum up the four terms
  w1 += shiftField(v1, mu, Backward());

  for(int site = 0; site < Nvol_; ++site)
    SetMat(field_strength, anti_hermite(mat(w1,site)), site);

  field_strength *= 0.25;
}

/////////////////////////////////////////////////////////////////////
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 1, \nu = 2\f$ */
void Dirac_Clover::isigma_12(FermionField& w,const FermionField& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site),-v.data[ff_.index_i(c,0,site)]);
      w.data.set(ff_.index_i(c,0,site), v.data[ff_.index_r(c,0,site)]);
      w.data.set(ff_.index_r(c,1,site), v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_i(c,1,site),-v.data[ff_.index_r(c,1,site)]);
      w.data.set(ff_.index_r(c,2,site),-v.data[ff_.index_i(c,2,site)]);
      w.data.set(ff_.index_i(c,2,site), v.data[ff_.index_r(c,2,site)]);
      w.data.set(ff_.index_r(c,3,site), v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_i(c,3,site),-v.data[ff_.index_r(c,3,site)]);
    }
  }
}
void Dirac_Clover::isigma_13(FermionField& w,const FermionField& v) const{
  isigma_31(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_14(FermionField& w,const FermionField& v) const{
  isigma_41(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_21(FermionField& w,const FermionField& v) const{
  isigma_12(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_23(FermionField& w,const FermionField& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site),-v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_i(c,0,site), v.data[ff_.index_r(c,1,site)]);
      w.data.set(ff_.index_r(c,1,site),-v.data[ff_.index_i(c,0,site)]);
      w.data.set(ff_.index_i(c,1,site), v.data[ff_.index_r(c,0,site)]);
      w.data.set(ff_.index_r(c,2,site),-v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_i(c,2,site), v.data[ff_.index_r(c,3,site)]);
      w.data.set(ff_.index_r(c,3,site),-v.data[ff_.index_i(c,2,site)]);
      w.data.set(ff_.index_i(c,3,site), v.data[ff_.index_r(c,2,site)]);
    }
  }
}
void Dirac_Clover::isigma_24(FermionField& w,const FermionField& v) const{
  isigma_42(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_31(FermionField& w,const FermionField& v) const{

  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site), v.data[ff_.index_r(c,1,site)]);
      w.data.set(ff_.index_i(c,0,site), v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_r(c,1,site),-v.data[ff_.index_r(c,0,site)]);
      w.data.set(ff_.index_i(c,1,site),-v.data[ff_.index_i(c,0,site)]);
      w.data.set(ff_.index_r(c,2,site), v.data[ff_.index_r(c,3,site)]);
      w.data.set(ff_.index_i(c,2,site), v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_r(c,3,site),-v.data[ff_.index_r(c,2,site)]);
      w.data.set(ff_.index_i(c,3,site),-v.data[ff_.index_i(c,2,site)]);
    }
  }
}
void Dirac_Clover::isigma_32(FermionField& w,const FermionField& v) const{
  isigma_23(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_34(FermionField& w,const FermionField& v) const{
  isigma_43(w,v);
  w *= -1.0;
}
void Dirac_Clover::isigma_41(FermionField& w,const FermionField& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site), v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_i(c,0,site),-v.data[ff_.index_r(c,3,site)]);
      w.data.set(ff_.index_r(c,1,site), v.data[ff_.index_i(c,2,site)]);
      w.data.set(ff_.index_i(c,1,site),-v.data[ff_.index_r(c,2,site)]);
      w.data.set(ff_.index_r(c,2,site), v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_i(c,2,site),-v.data[ff_.index_r(c,1,site)]);
      w.data.set(ff_.index_r(c,3,site), v.data[ff_.index_i(c,0,site)]);
      w.data.set(ff_.index_i(c,3,site),-v.data[ff_.index_r(c,0,site)]);
    }
  }
}
void Dirac_Clover::isigma_42(FermionField& w,const FermionField& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site),-v.data[ff_.index_r(c,3,site)]);
      w.data.set(ff_.index_i(c,0,site),-v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_r(c,1,site), v.data[ff_.index_r(c,2,site)]);
      w.data.set(ff_.index_i(c,1,site), v.data[ff_.index_i(c,2,site)]);
      w.data.set(ff_.index_r(c,2,site),-v.data[ff_.index_r(c,1,site)]);
      w.data.set(ff_.index_i(c,2,site),-v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_r(c,3,site), v.data[ff_.index_r(c,0,site)]);
      w.data.set(ff_.index_i(c,3,site), v.data[ff_.index_i(c,0,site)]);
    }
  }
}
void Dirac_Clover::isigma_43(FermionField& w,const FermionField& v) const{
  for(int site=0; site<Nvol_; ++site){
    for(int c=0; c<NC_; ++c){
      w.data.set(ff_.index_r(c,0,site), v.data[ff_.index_i(c,2,site)]);
      w.data.set(ff_.index_i(c,0,site),-v.data[ff_.index_r(c,2,site)]);
      w.data.set(ff_.index_r(c,1,site),-v.data[ff_.index_i(c,3,site)]);
      w.data.set(ff_.index_i(c,1,site), v.data[ff_.index_r(c,3,site)]);
      w.data.set(ff_.index_r(c,2,site), v.data[ff_.index_i(c,0,site)]);
      w.data.set(ff_.index_i(c,2,site),-v.data[ff_.index_r(c,0,site)]);
      w.data.set(ff_.index_r(c,3,site),-v.data[ff_.index_i(c,1,site)]);
      w.data.set(ff_.index_i(c,3,site), v.data[ff_.index_r(c,1,site)]);
    }
  }
}

const Field Dirac_Clover::gamma5(const Field& f) const{
  return Dw->gamma5(f);
}
const Field Dirac_Clover::mult(const Field& f) const{
  FermionField w, w2;

  w.data  = Dw->mult(f);
  mult_sw(w2,FermionField(f));
  w -= w2;
  return w.data;
}

const Field Dirac_Clover::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

//====================================================================
const Field Dirac_Clover::md_force(const Field& eta,const Field& zeta)const{
  //Wilson term
  Field force = Dw->md_force(eta,zeta);

 //just temporaries here (to be eliminated when all Field->FermionField)
  FermionField eta_F(eta), zeta_F(zeta);
 
  force += md_force_block(eta_F, zeta_F);
  force += md_force_block(zeta_F, eta_F);

  return force;
}
//====================================================================
const Field Dirac_Clover::md_force_block(const FermionField& eta,
					 const FermionField& zeta)const{
  using namespace SUNmatUtils;
  using namespace SUNvec_utils;
  using namespace FieldUtils;
  using namespace Mapping;

  //.................... Temporary variables declaration
  int Nd = CommonPrms::instance()->Nd();

  GaugeField force;
  GaugeField1D fce_tmp1, fce_tmp2;
  FermionField eta2, vleft, vright;
  FermionField shifted;
  GaugeField1D U_mu, U_nu;
  GaugeField1D Cmu_up, Cnu_up;
  Staples stpl_;
  SUNvec temp_v;
  //...................................................
  //really temporary
  GaugeField temp_u_(*u_);  
  
  //Clover term: 8 terms, see Matsufuru-san's note

  for(int mu=0; mu<NDIM_; ++mu){
    for(int nu=0; nu<NDIM_; ++nu){
      if(nu == mu) continue;

      (this->*isigma[NDIM_*mu+nu])(eta2,eta); //sigma_{mu,nu} mult

      U_nu = DirSlice(temp_u_,nu);
      U_mu = DirSlice(temp_u_,mu);

      //term 1 and 5 --------------------------------------
      Cmu_up  = stpl_.upper(temp_u_,mu,nu); //V_mu
      Cmu_up -= stpl_.lower(temp_u_,mu,nu); //V_+mu - V_-mu 
      for(int site = 0; site<Nvol_; ++site){
	for(int s=0; s<Nd; ++s){
	  temp_v = mat_dag(Cmu_up,site)*vec(eta2,s,site);
	  temp_v = mat(U_mu,site)*temp_v;
	  SetVec(vright, temp_v, s, site);
	} 
      }
      fce_tmp2 = field_oprod(zeta,vright); 

      //term 2 --------------------------------------
      Cnu_up = stpl_.upper(temp_u_, nu, mu);
      shifted = shiftField(eta2, nu, Forward());
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(Cnu_up,site))* vec(shifted,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      }
      shifted = shiftField(zeta, nu, Forward());
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = mat(U_nu,site) * vec(shifted,s,site);
	  SetVec(vleft, temp_v, s, site);
	}
      }// U_nu(x)*zeta(x+nu)
      fce_tmp2 += field_oprod(vleft,vright);

      //term 4 and 8 --------------------------------------
      shifted = shiftField(eta2, mu, Forward());
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = mat(U_mu,site) * vec(shifted,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      } // U_mu (x) * eta2(x+mu)   
      shifted = shiftField(zeta, mu, Forward());
      for(int site = 0; site < Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = mat(Cmu_up,site) * vec(shifted,s,site);
	  SetVec(vleft, temp_v, s, site);
	}
      } //  (V_+mu + V_-mu) * zeta(x+mu)
      fce_tmp2 += field_oprod(vleft,vright);
      
      //term 3 --------------------------------------
      shifted = shiftField(eta2, nu, Forward());
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  //temp_v = (mat(U_nu,site))* vec(shifted,s,site);
	  SetVec(vright, (mat(U_nu,site)* vec(shifted,s,site)), s, site);
	}
      }       
      shifted = shiftField(vright, mu, Forward()); //U_nu(x+mu)*eta2(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(U_mu,site))* vec(shifted,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      } //vright = U_mu(x)*U_nu(x+mu)*eta2(x+mu+nu)  
     
      shifted = shiftField(zeta, mu, Forward());
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(U_mu,site))* vec(shifted,s,site);
	  SetVec(vleft, temp_v, s, site);
	}
      }        
      shifted = shiftField(vleft, nu, Forward()); //U_mu(x+nu)*zeta(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(U_nu,site))* vec(shifted,s,site);
	  SetVec(vleft, temp_v, s, site);
	}
      } //vright = U_nu(x)*U_mu(x+nu)*zeta(x+mu+nu)  
      fce_tmp2 += field_oprod(vleft,vright);

      //term 6 --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat_dag(Cnu_up,site))* vec(eta2,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      }       
      shifted = shiftField(vright, nu, Backward()); // V_+nu (x-nu)*eta2(x-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat_dag(U_nu,site))* vec(zeta,s,site);
	  SetVec(vleft, temp_v, s, site);
	}
      }   
      fce_tmp2 -= field_oprod(shiftField(vleft, nu, Backward()), shifted);

      // term 7  --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat_dag(U_nu,site))* vec(eta2,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      }   
      shifted = shiftField(vright,nu,Backward());//shifted = Udag_nu(x-nu)*eta2(x-nu)
      shifted = shiftField(shifted,mu,Forward());//Udag_nu(x+mu-nu)*eta2(x+mu-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(U_mu,site))* vec(shifted,s,site);
	  SetVec(vright, temp_v, s, site);
	}
      } //vright = U_mu(x)*Udag_nu(x+mu-nu)*eta2(x+mu-nu)
     
      shifted = shiftField(zeta, mu, Forward()); //zeta(x+mu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (mat(U_mu,site))* vec(shifted,s,site);
	  temp_v = (mat_dag(U_nu,site))*temp_v;
	  SetVec(vleft, temp_v, s, site);
	}
      }       
      //shift = Udag_nu(x-nu)*U_mu(x-nu)*zeta(x+mu-nu)
      fce_tmp2 -= field_oprod(shiftField(vleft, nu, Backward()), vright);
      fce_tmp2 *= -Dw->getKappa()*csw_/8.0; 

      for(int site=0; site<Nvol_; ++site)
	AddMat(force,mat(fce_tmp2,site),site,mu); 
    }
  }
  return force.data;
}
//====================================================================
const vector<int> Dirac_Clover::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}
