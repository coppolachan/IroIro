/*! @file cloverTerm.cpp
 *  @brief implementation of the CloverTerm class
 *  Time-stamp: <2014-09-06 08:08:56 noaki>
 */
#include "cloverTerm.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include "include/messages_macros.hpp"

using namespace std;

void CloverTerm::mat_vec(double* wp,double* Fp,double* vp)const{
  for(int s=0; s<ND_; ++s){
    for(int c=0; c<NC_; ++c){
      double wr=0.0, wi=0.0;
      for(int c1=0; c1<NC_; ++c1){
	wr += Fp[re(c,c1)]*vp[fr(s,c1)] -Fp[im(c,c1)]*vp[fi(s,c1)];
	wi += Fp[re(c,c1)]*vp[fi(s,c1)] +Fp[im(c,c1)]*vp[fr(s,c1)];
      }
      wp[fr(s,c)] += wr;
      wp[fi(s,c)] += wi;
    }
  }
}

/*
void (CloverTerm::*CloverTerm::isigma[])(FermionField&,
					 const FermionField&)const
= {&CloverTerm::isigma_diag,
   &CloverTerm::isigma_12,
   &CloverTerm::isigma_13,
   &CloverTerm::isigma_14,
   &CloverTerm::isigma_21,
   &CloverTerm::isigma_diag,
   &CloverTerm::isigma_23,
   &CloverTerm::isigma_24,
   &CloverTerm::isigma_31,
   &CloverTerm::isigma_32,
   &CloverTerm::isigma_diag,
   &CloverTerm::isigma_34,
   &CloverTerm::isigma_41,
   &CloverTerm::isigma_42,
   &CloverTerm::isigma_43,
   &CloverTerm::isigma_diag,};

void CloverTerm::isigma_diag(Field&,const Field&) const{
  CCIO::cout << "Illegal value of (mu,nu) in CloverTerm::isigma.\n";
  abort();
}
*/

void CloverTerm::mult_sw(Field& w,const Field& v)const{

  for(int site=0; site<Nvol_; ++site){
    double* vp = const_cast<Field&>(v).getaddr(ff_.index(0,site));
    double* wp = w.getaddr(ff_.index(0,site));
    vector<double> wt(ff_.Nin());

    dm_.isigma23core(&wt[0],vp);
    mat_vec(wp,Bx_.data.getaddr(gf_.index(0,site)),&wt[0]);
    dm_.isigma31core(&wt[0],vp);
    mat_vec(wp,By_.data.getaddr(gf_.index(0,site)),&wt[0]);
    dm_.isigma12core(&wt[0],vp);
    mat_vec(wp,Bz_.data.getaddr(gf_.index(0,site)),&wt[0]);
    dm_.isigma41core(&wt[0],vp);
    mat_vec(wp,Ex_.data.getaddr(gf_.index(0,site)),&wt[0]);
    dm_.isigma42core(&wt[0],vp);
    mat_vec(wp,Ey_.data.getaddr(gf_.index(0,site)),&wt[0]);
    dm_.isigma43core(&wt[0],vp);
    mat_vec(wp,Ez_.data.getaddr(gf_.index(0,site)),&wt[0]);
  }
  w *= csw_;
}

void CloverTerm::set_csw() {
  _Message(DEBUG_VERB_LEVEL, "[DiracClover] Setting Field Strenght\n");
  Staples  stpl;
  Bx_= stpl.fieldStrength(GaugeField(*u_),1,2);
  By_= stpl.fieldStrength(GaugeField(*u_),2,0);
  Bz_= stpl.fieldStrength(GaugeField(*u_),0,1);
  Ex_= stpl.fieldStrength(GaugeField(*u_),3,0);
  Ey_= stpl.fieldStrength(GaugeField(*u_),3,1);
  Ez_= stpl.fieldStrength(GaugeField(*u_),3,2);
}

const Field CloverTerm::mult(const Field& f)const{
  Field w(f.size());
  mult_sw(w,f);
  return w;
}

/*
//====================================================================
const Field CloverTerm::md_force(const Field& eta,const Field& zeta)const{
  //Wilson term
  Field force = Dw_->md_force(eta,zeta);

 //just temporaries here (to be eliminated when all Field->FermionField)
  FermionField eta_F(eta), zeta_F(zeta);
 
  force += md_force_block(eta_F, zeta_F);
  force += md_force_block(zeta_F, eta_F);

  return force;
}
//====================================================================
const Field CloverTerm::md_force_block(const FermionField& eta,
					 const FermionField& zeta)const{
  using namespace SUNmatUtils;
  using namespace SUNvecUtils;
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
  GaugeField temp_u(*u_);  
  
  //Clover term: 8 terms, see Matsufuru-san's note

  for(int mu=0; mu<NDIM_; ++mu){
    for(int nu=0; nu<NDIM_; ++nu){
      if(nu == mu) continue;

      (this->*isigma[NDIM_*mu+nu])(eta2,eta); //sigma_{mu,nu} mult

      U_nu = DirSlice(temp_u,nu);
      U_mu = DirSlice(temp_u,mu);

      //term 1 and 5 --------------------------------------
      Cmu_up  = stpl_.upper(temp_u,mu,nu); //V_mu
      Cmu_up -= stpl_.lower(temp_u,mu,nu); //V_+mu - V_-mu 
      for(int site = 0; site<Nvol_; ++site){
	for(int s=0; s<Nd; ++s){
	  temp_v = mat_dag(Cmu_up,site)*vec(eta2,s,site);
	  temp_v = mat(U_mu,site)*temp_v;
	  SetVec(vright, temp_v, s, site);
	} 
      }
      fce_tmp2 = field_oprod(zeta,vright); 

      //term 2 --------------------------------------
      Cnu_up = stpl_.upper(temp_u, nu, mu);
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
      fce_tmp2 *= -kpp_*csw_/8.0; 

      for(int site=0; site<Nvol_; ++site)
	AddMat(force,mat(fce_tmp2,site),site,mu); 
    }
  }
  return force.data;
}
*/
