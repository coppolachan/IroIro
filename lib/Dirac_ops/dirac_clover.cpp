/*! 
  @file dirac_clover.cpp

  @brief Dirac Clover operator class. Functions definition
*/

#include "dirac_clover.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Communicator/comm_io.hpp"
#include "Main/Geometry/mapper.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/messages_macros.hpp"

using namespace std;

//=================================================================

//void (Dirac_Clover::*Dirac_Clover::mult_isigma[])(Field&,const Field&)const 
//= {&Dirac_Clover::i_sigma};
  
void Dirac_Clover::mult_isigma(FermionField& v, const FermionField& w,int mu,int nu) const {
  
  if(mu==nu) {
    CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_isigma.\n";
    abort();
  }

  if(mu == 0){
    switch(nu) {
    case (1): 
      mult_isigma12(v,w);
      break;
    case (2):
      mult_isigma31(v,w);
      v *= -1.0;
      break;
    case (3):
      mult_isigma41(v,w);
      v *= -1.0;
      break;
    default:
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_isigma.\n";
      abort();
    }
  }

  if(mu == 1) {
    switch(nu) {
    case(0):
      mult_isigma12(v,w);
      v *= -1.0;
      break;
    case(2):
      mult_isigma23(v,w);
      break;
    case(3):
      mult_isigma42(v,w);
      v *= -1.0; 
      break;
    default:
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_isigma.\n";
      abort();
    }
  }
  
  if(mu == 2) {
    switch(nu) {
    case (0):
      mult_isigma31(v,w);
      break;
    case(1):
      mult_isigma23(v,w);
      v *= -1.0;
      break;
    case(3):
      mult_isigma43(v,w);
      v *= -1.0;
      break;
    default:
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_isigma.\n";
      abort();
    }
  }
  

  if(mu == 3){
    switch(nu) {
    case(0):
      mult_isigma41(v,w);
      break;
    case(1):
      mult_isigma42(v,w); 
      break;
    case(2):
      mult_isigma43(v,w);
      break;
    default:
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_isigma.\n";
      abort();
    }
  }
  
  
}



//====================================================================
void Dirac_Clover::mult_csw(FermionField& v_out, const FermionField& w) const {
  using namespace SUNmatUtils;
  using namespace SUNvec_utils;
  using namespace FieldUtils;

  //Temporary variable declarations
  SUNvec v1;
  FermionField wt;

  //Actual calculation - START
  mult_isigma23(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      SetVector(v_out, (matrix(d_Bx,site) * vect(wt,s,site)), s, site);
    }
  }
   
  mult_isigma31(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      AddVector(v_out, (matrix(d_By,site) * vect(wt,s,site)), s, site);
    }
  }

  mult_isigma12(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      AddVector(v_out, (matrix(d_Bz,site) * vect(wt,s,site)), s, site);
    }
  }

  mult_isigma41(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      AddVector(v_out, (matrix(d_Ex,site) * vect(wt,s,site)), s, site);
    }
  }

  mult_isigma42(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      AddVector(v_out, (matrix(d_Ey,site) * vect(wt,s,site)), s, site);
    }
  }
  

  mult_isigma43(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < NDIM_; ++s){
      v1 = matrix(d_Ez,site) * vect(wt,s,site);
      AddVector(v_out, v1, s, site);
    }
  }
  
  v_out *= Dw->getKappa() * csw_;
}


//====================================================================
void Dirac_Clover::set_csw() {

  _Message(DEBUG_VERB_LEVEL, "[DiracClover] Setting Field Strenght\n");

  set_fieldstrength(d_Bx, 1, 2);
  set_fieldstrength(d_By, 2, 0);
  set_fieldstrength(d_Bz, 0, 1);
  set_fieldstrength(d_Ex, 3, 0);
  set_fieldstrength(d_Ey, 3, 1);
  set_fieldstrength(d_Ez, 3, 2);

}

//====================================================================
/*! @brief Calculates the term \f$F_{\mu.\nu}\f$  */
void Dirac_Clover::set_fieldstrength(GaugeField1D& field_strength,
				     const int mu, const int nu){
  using namespace SUNmatUtils;
  using namespace FieldUtils;
  using namespace MapsEnv;

  //.................. Temporary variables declaration
  GaugeField1D Cup, Cdn;
  GaugeField1D U_mu;  /*< @brief \f$U_\mu(x)\f$ */
  GaugeField1D w1, w2, v1, v2;
  Staples stpl_;

  //really temporary
  GaugeField temp_u_(*u_);
  //.......................................
  U_mu = DirSlice(temp_u_, mu);

  Cup = stpl_.upper(temp_u_,mu,nu); // Upper staple V_+mu
  Cdn = stpl_.lower(temp_u_,mu,nu); // Lower staple V_-mu

  for(int site = 0; site < Nvol_; ++site){
    SetMatrix(w1, matrix(U_mu,site)    * matrix_dag(Cup,site), site);// U_mu(x)*(V_+mu)^dag
    SetMatrix(w2, matrix(U_mu,site)    * matrix_dag(Cdn,site), site);// U_mu(x)*(V_-mu)^dag
    SetMatrix(v1, matrix_dag(Cup,site) * matrix(U_mu,site)   , site);// (V_+mu)^dag*U_mu(x)
    SetMatrix(v2, matrix_dag(Cdn,site) * matrix(U_mu,site)   , site);// (V_-mu)^dag*U_mu(x)
  }

  w1 -= w2;

  //    +--<--+ 
  //    |     | w1
  //    |     |
  // (x)+-->--+
  //    |     |
  //    |     | w2
  //    +--<--+  

  v1 -= v2; 
 
  //    +--<--+ 
  //    |     | v1
  //    |     |
  //    +-->--+(x)
  //    |     |
  //    |     | v2
  //    +--<--+  

  //Sum up the four terms
  w1 += shift(v1, mu, Backward);

  for(int site = 0; site < Nvol_; ++site){
    SetMatrix(field_strength, anti_hermite(matrix(w1,site)), site);
  }

  field_strength *= 0.25;

}
//====================================================================
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 2, \nu = 3\f$ */
void Dirac_Clover::mult_isigma23(FermionField& v, const FermionField& w) const {
  // v = \sigma_23 * w

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site), -w.data[ff_->index_i(cc,1,site)]);
      v.data.set(ff_->index_i(cc,0,site),  w.data[ff_->index_r(cc,1,site)]);
      
      v.data.set(ff_->index_r(cc,1,site), -w.data[ff_->index_i(cc,0,site)]);
      v.data.set(ff_->index_i(cc,1,site),  w.data[ff_->index_r(cc,0,site)]);
      
      v.data.set(ff_->index_r(cc,2,site), -w.data[ff_->index_i(cc,3,site)]);
      v.data.set(ff_->index_i(cc,2,site),  w.data[ff_->index_r(cc,3,site)]);
      
      v.data.set(ff_->index_r(cc,3,site), -w.data[ff_->index_i(cc,2,site)]);
      v.data.set(ff_->index_i(cc,3,site),  w.data[ff_->index_r(cc,2,site)]);
      
    }
  }
}
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 3, \nu = 1\f$ */
//====================================================================
void Dirac_Clover::mult_isigma31(FermionField& v, const FermionField& w) const {

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site),  w.data[ff_->index_r(cc,1,site)]);
      v.data.set(ff_->index_i(cc,0,site),  w.data[ff_->index_i(cc,1,site)]);
      
      v.data.set(ff_->index_r(cc,1,site), -w.data[ff_->index_r(cc,0,site)]);
      v.data.set(ff_->index_i(cc,1,site), -w.data[ff_->index_i(cc,0,site)]);
      
      v.data.set(ff_->index_r(cc,2,site),  w.data[ff_->index_r(cc,3,site)]);
      v.data.set(ff_->index_i(cc,2,site),  w.data[ff_->index_i(cc,3,site)]);
      
      v.data.set(ff_->index_r(cc,3,site), -w.data[ff_->index_r(cc,2,site)]);
      v.data.set(ff_->index_i(cc,3,site), -w.data[ff_->index_i(cc,2,site)]);
      
    }
  }
}
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 1, \nu = 2\f$ */
//====================================================================
void Dirac_Clover::mult_isigma12(FermionField& v, const FermionField& w) const {

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site), -w.data[ff_->index_i(cc,0,site)]);
      v.data.set(ff_->index_i(cc,0,site),  w.data[ff_->index_r(cc,0,site)]);
      
      v.data.set(ff_->index_r(cc,1,site),  w.data[ff_->index_i(cc,1,site)]);
      v.data.set(ff_->index_i(cc,1,site), -w.data[ff_->index_r(cc,1,site)]);
      
      v.data.set(ff_->index_r(cc,2,site), -w.data[ff_->index_i(cc,2,site)]);
      v.data.set(ff_->index_i(cc,2,site),  w.data[ff_->index_r(cc,2,site)]);
      
      v.data.set(ff_->index_r(cc,3,site),  w.data[ff_->index_i(cc,3,site)]);
      v.data.set(ff_->index_i(cc,3,site), -w.data[ff_->index_r(cc,3,site)]);
      
    }
  }
}
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 4, \nu = 1\f$ */
//====================================================================
void Dirac_Clover::mult_isigma41(FermionField& v, const FermionField& w) const {
  
  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site),  w.data[ff_->index_i(cc,3,site)]);
      v.data.set(ff_->index_i(cc,0,site), -w.data[ff_->index_r(cc,3,site)]);
      
      v.data.set(ff_->index_r(cc,1,site),  w.data[ff_->index_i(cc,2,site)]);
      v.data.set(ff_->index_i(cc,1,site), -w.data[ff_->index_r(cc,2,site)]);
      
      v.data.set(ff_->index_r(cc,2,site),  w.data[ff_->index_i(cc,1,site)]);
      v.data.set(ff_->index_i(cc,2,site), -w.data[ff_->index_r(cc,1,site)]);
      
      v.data.set(ff_->index_r(cc,3,site),  w.data[ff_->index_i(cc,0,site)]);
      v.data.set(ff_->index_i(cc,3,site), -w.data[ff_->index_r(cc,0,site)]);
      
    }
  }
}
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 4, \nu = 2\f$ */
//====================================================================
void Dirac_Clover::mult_isigma42(FermionField& v, const FermionField& w) const {

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site), -w.data[ff_->index_r(cc,3,site)]);
      v.data.set(ff_->index_i(cc,0,site), -w.data[ff_->index_i(cc,3,site)]);
      
      v.data.set(ff_->index_r(cc,1,site),  w.data[ff_->index_r(cc,2,site)]);
      v.data.set(ff_->index_i(cc,1,site),  w.data[ff_->index_i(cc,2,site)]);

      v.data.set(ff_->index_r(cc,2,site), -w.data[ff_->index_r(cc,1,site)]);
      v.data.set(ff_->index_i(cc,2,site), -w.data[ff_->index_i(cc,1,site)]);
      
      v.data.set(ff_->index_r(cc,3,site),  w.data[ff_->index_r(cc,0,site)]);
      v.data.set(ff_->index_i(cc,3,site),  w.data[ff_->index_i(cc,0,site)]);
      
    }
  }
}
/*! @brief Calculates the product \f$\sigma_{\mu,\nu} v\f$, 
  \f$\mu = 4, \nu = 3\f$ */
//====================================================================
void Dirac_Clover::mult_isigma43(FermionField& v, const FermionField& w) const {

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < NC_; ++cc){
      
      v.data.set(ff_->index_r(cc,0,site),  w.data[ff_->index_i(cc,2,site)]);
      v.data.set(ff_->index_i(cc,0,site), -w.data[ff_->index_r(cc,2,site)]);
      
      v.data.set(ff_->index_r(cc,1,site), -w.data[ff_->index_i(cc,3,site)]);
      v.data.set(ff_->index_i(cc,1,site),  w.data[ff_->index_r(cc,3,site)]);
      
      v.data.set(ff_->index_r(cc,2,site),  w.data[ff_->index_i(cc,0,site)]);
      v.data.set(ff_->index_i(cc,2,site), -w.data[ff_->index_r(cc,0,site)]);
      
      v.data.set(ff_->index_r(cc,3,site), -w.data[ff_->index_i(cc,1,site)]);
      v.data.set(ff_->index_i(cc,3,site),  w.data[ff_->index_r(cc,1,site)]);
      
    }
  }
}
//====================================================================
const Field Dirac_Clover::gamma5(const Field& f) const{
  return Dw->gamma5(f);
}
//====================================================================
const Field Dirac_Clover::mult(const Field& f) const{
  FermionField w, w2;

  w.data  = Dw->mult(f);
  mult_csw(w2,FermionField(f));
  w -= w2;
  return w.data;
}
//====================================================================
const Field Dirac_Clover::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}
//====================================================================
//temporarily here

/*! @brief Calculates the external product of two vectors

  \f[(A^\dagger \wedge B)_{ab} = A^*_a B_b \f]
 */
void Dirac_Clover::external_prod(GaugeField1D& res, 
				 const FermionField& A, 
				 const FermionField& B) const {
  using namespace FieldUtils;
  assert(A.size() == B.size());

  // .................. Temporary variables declaration
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;
  double fre, fim;
  size_t ra, rb, ia, ib; // for index calculations
  //...................................................

  for(int site = 0; site < Nvol_; ++site){
    f = 0.0;
    for(int a = 0; a < NC_; ++a){
      for(int b = 0; b < NC_; ++b){
	fre = fim = 0.0;
	for(int s=0; s<Nd; ++s){
	  //indexes
	  ra =A.format.index_r(a,s,site);
	  ia =A.format.index_i(a,s,site);
	  
	  rb =B.format.index_r(b,s,site);
	  ib =B.format.index_i(b,s,site);
	  
	  fre += A.data[rb]*B.data[ra] + A.data[ib]*B.data[ia];
	  fim += A.data[rb]*B.data[ia] - A.data[ib]*B.data[ra];
	}
	f.set(a,b,fre,fim);
      }
    }
    // Store matrix in res field
    SetMatrix(res, f, gsite(site));
  }
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
  using namespace MapsEnv;

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

  for(int mu = 0; mu < NDIM_; ++mu){
    for(int nu = 0; nu < NDIM_; ++nu){
      if(nu == mu) continue;

      mult_isigma(eta2, eta, mu, nu); //sigma_{mu,nu} mult

      U_nu = DirSlice(temp_u_, nu);
      U_mu = DirSlice(temp_u_, mu);

      //term 1 and 5 --------------------------------------
      Cmu_up  = stpl_.upper(temp_u_, mu,nu); //V_mu
      Cmu_up -= stpl_.lower(temp_u_, mu,nu); //V_+mu - V_-mu 
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = matrix_dag(Cmu_up,site) * vect(eta2,s,site);
	  temp_v = matrix(U_mu,site) * temp_v;
	  SetVector(vright, temp_v, s, site);
	} 
      }
      external_prod(fce_tmp1, zeta, vright);
      fce_tmp2 = fce_tmp1;
 
  
 
      //term 2 --------------------------------------
      Cnu_up = stpl_.upper(temp_u_, nu, mu);
      shifted = shift(eta2, nu, Forward);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(Cnu_up,site))* vect(shifted,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      }
      shifted = shift(zeta, nu, Forward);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = matrix(U_nu,site) * vect(shifted,s,site);
	  SetVector(vleft, temp_v, s, site);
	}
      }// U_nu(x)*zeta(x+nu)
      external_prod(fce_tmp1, vleft, vright);
      fce_tmp2 += fce_tmp1;

   

      //term 4 and 8 --------------------------------------
      shifted = shift(eta2, mu, Forward);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = matrix(U_mu,site) * vect(shifted,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      } // U_mu (x) * eta2(x+mu)   
      shifted = shift(zeta, mu, Forward);
      for(int site = 0; site < Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = matrix(Cmu_up,site) * vect(shifted,s,site);
	  SetVector(vleft, temp_v, s, site);
	}
      } //  (V_+mu + V_-mu) * zeta(x+mu)
      external_prod(fce_tmp1, vleft, vright);
      fce_tmp2 += fce_tmp1;

  

      //term 3 --------------------------------------
      shifted = shift(eta2, nu, Forward);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  //temp_v = (matrix(U_nu,site))* vect(shifted,s,site);
	  SetVector(vright, (matrix(U_nu,site)* vect(shifted,s,site)), s, site);
	}
      }       
      shifted = shift(vright, mu, Forward); //U_nu(x+mu)*eta2(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(U_mu,site))* vect(shifted,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      } //vright = U_mu(x)*U_nu(x+mu)*eta2(x+mu+nu)  
     
      shifted = shift(zeta, mu, Forward);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(U_mu,site))* vect(shifted,s,site);
	  SetVector(vleft, temp_v, s, site);
	}
      }        
      shifted = shift(vleft, nu, Forward); //U_mu(x+nu)*zeta(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(U_nu,site))* vect(shifted,s,site);
	  SetVector(vleft, temp_v, s, site);
	}
      } //vright = U_nu(x)*U_mu(x+nu)*zeta(x+mu+nu)  
      external_prod(fce_tmp1, vleft, vright);
      fce_tmp2 += fce_tmp1;

 
      //term 6 --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix_dag(Cnu_up,site))* vect(eta2,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      }       
      shifted = shift(vright, nu, Backward); // V_+nu (x-nu)*eta2(x-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix_dag(U_nu,site))* vect(zeta,s,site);
	  SetVector(vleft, temp_v, s, site);
	}
      }   
      external_prod(fce_tmp1, shift(vleft, nu, Backward), shifted);
      fce_tmp2 -= fce_tmp1;   

      // term 7  --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix_dag(U_nu,site))* vect(eta2,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      }   
      shifted = shift(vright, nu, Backward); // shifted = Udag_nu(x-nu)*eta2(x-nu)
      shifted = shift(shifted, mu, Forward); // Udag_nu(x+mu-nu)*eta2(x+mu-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(U_mu,site))* vect(shifted,s,site);
	  SetVector(vright, temp_v, s, site);
	}
      } //vright = U_mu(x)*Udag_nu(x+mu-nu)*eta2(x+mu-nu)
     
      shifted = shift(zeta, mu, Forward); //zeta(x+mu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  temp_v = (matrix(U_mu,site))* vect(shifted,s,site);
	  temp_v = (matrix_dag(U_nu,site))*temp_v;
	  SetVector(vleft, temp_v, s, site);
	}
      }       
      //shift = Udag_nu(x-nu)*U_mu(x-nu)*zeta(x+mu-nu)
      external_prod(fce_tmp1, shift(vleft, nu, Backward), vright);
      fce_tmp2 -= fce_tmp1;   
     
      fce_tmp2 *= - Dw->getKappa() * csw_ / 8.0; 

      for(int site = 0; site<Nvol_; ++site)
	AddMatrix(force, matrix(fce_tmp2, site), site, mu); 
      


    }
  }

  return force.data;
}
//====================================================================
const vector<int> Dirac_Clover::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}
