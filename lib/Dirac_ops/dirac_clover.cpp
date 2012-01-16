/*! 
  @file dirac_clover.cpp

  @brief Dirac Clover operator class. Functions definition
*/

#include "dirac_clover.hpp"
#include "Communicator/comm_io.hpp"
#include "Tools/sunMatUtils.hpp"

using namespace std;

//======== Auxiliary routines

const std::valarray<double> Dirac_Clover::anti_herm(const SUNmat& m){
  
  std::valarray<double> va(m.getva());
  for(int a=0; a<NC_; ++a){
    va[2*(NC_*a+a)  ]= 0.0;
    va[2*(NC_*a+a)+1]= 0.0;
    
    for(int b=0; b<a; ++b){
      int ab = 2*(NC_*a+b);
      int ba = 2*(NC_*b+a);
      
      va[ab]-= va[ba];
      va[ab]/= 2.0;
      va[ba] =-va[ab];
      
      va[ab+1]+= va[ba+1];
      va[ab+1]/= 2.0;
      va[ba+1] = va[ab+1];
    }
  }
  return va;
}

//=================================================================

void Dirac_Clover::mult_isigma(Field& v, const Field& w,
			       const int mu, const int nu) const {
  
  if(mu==0){
    if(nu==1){
      mult_isigma12(v,w);
    }else if(nu==2){
      mult_isigma31(v,w);
      v *= -1.0;
    }else if(nu==3){
      mult_isigma41(v,w);
      v *= -1.0;
    }else{
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_sg.\n";
      abort();
    }
  }else if(mu==1){
    if(nu==0){
      mult_isigma12(v,w);
      v *= -1.0;
    }else if(nu==2){
      mult_isigma23(v,w);
    }else if(nu==3){
      mult_isigma42(v,w);
      v *= -1.0;
    }else{
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_sg.\n";
      abort();
    }
  }else if(mu==2){
    if(nu==0){
      mult_isigma31(v,w);
    }else if(nu==1){
      mult_isigma23(v,w);
      v *= -1.0;
    }else if(nu==3){
      mult_isigma43(v,w);
      v *= -1.0;
    }else{
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_sg.\n";
      abort();
    }
  }else if(mu==3){
    if(nu==0){
      mult_isigma41(v,w);
    }else if(nu==1){
      mult_isigma42(v,w);
    }else if(nu==2){
      mult_isigma43(v,w);
    }else{
      CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_sg.\n";
      abort();
    }
  }else{
    CCIO::cout << "Illegal value of (mu,nu) in Dirac_Clover::mult_sg.\n";
    abort();
  }

}



//====================================================================
void Dirac_Clover::mult_csw(Field& v_out, const Field& w) const {
  using namespace SUNmat_utils;
  using namespace SUNvec_utils;
  //  typedef std::valarray<double> vect;
  //vect v1(2*Nc);
  SUNvec v1;
  Field wt(fsize_);

  mult_isigma23(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_Bx,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site), v1.getva());
    }
  }

  mult_isigma31(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_By,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site), v1.getva());
    }
  }

  mult_isigma12(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_Bz,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site),v1.getva());
    }
  }
  
  mult_isigma41(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_Ex,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site), v1.getva());
    }
  }


  mult_isigma42(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_Ey,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site),v1.getva());
    }
  }
  

  mult_isigma43(wt,w);
  for(int site = 0; site < Nvol_; ++site){
    for(int s = 0; s < Ndim_; ++s){
      v1 = u(d_Ez,site) * v(wt,s,site);
      v_out.add(ff_->cslice(s,site), v1.getva());
    }
  }
 

  v_out *= Dw->getKappa() * csw_;
}

//====================================================================
void Dirac_Clover::set_csw() {

  set_fieldstrength(d_Bx, 1, 2);
  set_fieldstrength(d_By, 2, 0);
  set_fieldstrength(d_Bz, 0, 1);
  set_fieldstrength(d_Ex, 3, 0);
  set_fieldstrength(d_Ey, 3, 1);
  set_fieldstrength(d_Ez, 3, 2);

}

//====================================================================
void Dirac_Clover::set_fieldstrength(GaugeField1D& field_strength,
				     const int mu, const int nu){
  using namespace SUNmat_utils;

  GaugeField1D Cup, Cdn;
  GaugeField1D U_mu((*u_)[gf_->dir_slice(mu)]);  //U_mu
  GaugeField1D w, v1, v2;

  Cup.U = stpl_->upper(*u_,mu,nu);
  Cdn.U = stpl_->lower(*u_,mu,nu);

  for(int site = 0; site < Nvol_; ++site){
    w.set_matrix(site,u(U_mu,site)*u_dag(Cup,site));
    v2.set_matrix(site,u(U_mu,site)*u_dag(Cdn,site));
  }

  w.U -=v2.U;

  for(int site = 0; site < Nvol_; ++site){
    v1.set_matrix(site,u_dag(Cup,site)*u(U_mu,site));
    v2.set_matrix(site,u_dag(Cdn,site)*u(U_mu,site));
  }

  v1.U -= v2.U;  

  sf_up_[mu]->setf(v1.U);

  w.U += sf_up_[mu]->getva();

  for(int site = 0; site < Nvol_; ++site){
    field_strength.U.set(field_strength.Format.cslice(0,site,0),
			 anti_herm(u(w,site)));
  }

  field_strength.U *= 0.25;

}

//====================================================================
void Dirac_Clover::mult_isigma23(Field& v, const Field& w) const {
  // v = \sigma_23 * w

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site), -w[ff_->index_i(cc,1,site)]);
      v.set(ff_->index_i(cc,0,site),  w[ff_->index_r(cc,1,site)]);
      
      v.set(ff_->index_r(cc,1,site), -w[ff_->index_i(cc,0,site)]);
      v.set(ff_->index_i(cc,1,site),  w[ff_->index_r(cc,0,site)]);
      
      v.set(ff_->index_r(cc,2,site), -w[ff_->index_i(cc,3,site)]);
      v.set(ff_->index_i(cc,2,site),  w[ff_->index_r(cc,3,site)]);
      
      v.set(ff_->index_r(cc,3,site), -w[ff_->index_i(cc,2,site)]);
      v.set(ff_->index_i(cc,3,site),  w[ff_->index_r(cc,2,site)]);
      
    }
  }
  

}

//====================================================================
void Dirac_Clover::mult_isigma31(Field& v, const Field& w) const {

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site),  w[ff_->index_r(cc,1,site)]);
      v.set(ff_->index_i(cc,0,site),  w[ff_->index_i(cc,1,site)]);
      
      v.set(ff_->index_r(cc,1,site), -w[ff_->index_r(cc,0,site)]);
      v.set(ff_->index_i(cc,1,site), -w[ff_->index_i(cc,0,site)]);
      
      v.set(ff_->index_r(cc,2,site),  w[ff_->index_r(cc,3,site)]);
      v.set(ff_->index_i(cc,2,site),  w[ff_->index_i(cc,3,site)]);
      
      v.set(ff_->index_r(cc,3,site), -w[ff_->index_r(cc,2,site)]);
      v.set(ff_->index_i(cc,3,site), -w[ff_->index_i(cc,2,site)]);
      
    }
  }
  
  
}

//====================================================================
void Dirac_Clover::mult_isigma12(Field& v, const Field& w) const {

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site), -w[ff_->index_i(cc,0,site)]);
      v.set(ff_->index_i(cc,0,site),  w[ff_->index_r(cc,0,site)]);
      
      v.set(ff_->index_r(cc,1,site),  w[ff_->index_i(cc,1,site)]);
      v.set(ff_->index_i(cc,1,site), -w[ff_->index_r(cc,1,site)]);
      
      v.set(ff_->index_r(cc,2,site), -w[ff_->index_i(cc,2,site)]);
      v.set(ff_->index_i(cc,2,site),  w[ff_->index_r(cc,2,site)]);
      
      v.set(ff_->index_r(cc,3,site),  w[ff_->index_i(cc,3,site)]);
      v.set(ff_->index_i(cc,3,site), -w[ff_->index_r(cc,3,site)]);
      
    }
  }
  

}

//====================================================================
void Dirac_Clover::mult_isigma41(Field& v, const Field& w) const {

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site),  w[ff_->index_i(cc,3,site)]);
      v.set(ff_->index_i(cc,0,site), -w[ff_->index_r(cc,3,site)]);
      
      v.set(ff_->index_r(cc,1,site),  w[ff_->index_i(cc,2,site)]);
      v.set(ff_->index_i(cc,1,site), -w[ff_->index_r(cc,2,site)]);
      
      v.set(ff_->index_r(cc,2,site),  w[ff_->index_i(cc,1,site)]);
      v.set(ff_->index_i(cc,2,site), -w[ff_->index_r(cc,1,site)]);
      
      v.set(ff_->index_r(cc,3,site),  w[ff_->index_i(cc,0,site)]);
      v.set(ff_->index_i(cc,3,site), -w[ff_->index_r(cc,0,site)]);
      
    }
  }


}

//====================================================================
void Dirac_Clover::mult_isigma42(Field& v, const Field& w) const {

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site), -w[ff_->index_r(cc,3,site)]);
      v.set(ff_->index_i(cc,0,site), -w[ff_->index_i(cc,3,site)]);
      
      v.set(ff_->index_r(cc,1,site),  w[ff_->index_r(cc,2,site)]);
      v.set(ff_->index_i(cc,1,site),  w[ff_->index_i(cc,2,site)]);

      v.set(ff_->index_r(cc,2,site), -w[ff_->index_r(cc,1,site)]);
      v.set(ff_->index_i(cc,2,site), -w[ff_->index_i(cc,1,site)]);
      
      v.set(ff_->index_r(cc,3,site),  w[ff_->index_r(cc,0,site)]);
      v.set(ff_->index_i(cc,3,site),  w[ff_->index_i(cc,0,site)]);
      
    }
  }


}

//====================================================================
void Dirac_Clover::mult_isigma43(Field& v, const Field& w) const {

  int Nc = CommonPrms::Nc();

  for(int site = 0; site < Nvol_; ++site){
    for(int cc = 0; cc < Nc; ++cc){
      
      v.set(ff_->index_r(cc,0,site),  w[ff_->index_i(cc,2,site)]);
      v.set(ff_->index_i(cc,0,site), -w[ff_->index_r(cc,2,site)]);
      
      v.set(ff_->index_r(cc,1,site), -w[ff_->index_i(cc,3,site)]);
      v.set(ff_->index_i(cc,1,site),  w[ff_->index_r(cc,3,site)]);
      
      v.set(ff_->index_r(cc,2,site),  w[ff_->index_i(cc,0,site)]);
      v.set(ff_->index_i(cc,2,site), -w[ff_->index_r(cc,0,site)]);
      
      v.set(ff_->index_r(cc,3,site), -w[ff_->index_i(cc,1,site)]);
      v.set(ff_->index_i(cc,3,site),  w[ff_->index_r(cc,1,site)]);
      
    }
  }
  

}


const Field Dirac_Clover::gamma5(const Field& f) const{
  return Dw->gamma5(f);
}



const Field Dirac_Clover::mult(const Field& f) const{
  Field w (ff_->size());
  Field w2(ff_->size());

  w  = Dw->mult(f);
  mult_csw(w2,f);
  w -= w2;

  return w;
}

const Field Dirac_Clover::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

//temporary here
void Dirac_Clover::external_prod(Field& res, 
				 const Field& A, 
				 const Field& B) const {
  //external A^dag x B
  int Nd = CommonPrms::instance()->Nd();
  SUNmat f;

  for(int site=0; site< Nvol_; ++site){
    f=0.0;
    for(int a=0; a<NC_; ++a){
      for(int b=0; b<NC_; ++b){
	double fre = 0.0;
	double fim = 0.0;
	for(int s=0; s<Nd; ++s){
	  //indexes
	  size_t ra =ff_->index_r(a,s,site);
	  size_t ia =ff_->index_i(a,s,site);
	  
	  size_t rb =ff_->index_r(b,s,site);
	  size_t ib =ff_->index_i(b,s,site);
	  
	  fre += A[rb]*B[ra] + A[ib]*B[ia];
	  fim += A[rb]*B[ia] - A[ib]*B[ra];
	}
	f.set(a,b,fre,fim);
      }
    }
    res.set(gf_->cslice(0,gsite(site),0),f.getva());
  }
}

const Field Dirac_Clover::md_force(const Field& eta,const Field& zeta)const{
  using namespace SUNmat_utils;
  using namespace SUNvec_utils;
  int Nd = CommonPrms::instance()->Nd();

  GaugeField force;
  GaugeField1D fce_tmp1, fce_tmp2;
  Field eta2(fsize_), vleft(fsize_), vright(fsize_);
  Field shifted(fsize_);
  GaugeField1D U_mu, U_nu;
  GaugeField1D Cmu_up, Cnu_up;

  SUNvec vect;
  //valarray<double> vect(2*NC_);
  
  //Wilson term
  force.U = Dw->md_force(eta,zeta);

  //Clover term: 8 terms, see Matsufuru-san's note

  for(int mu = 0; mu < Ndim_; ++mu){
    for(int nu = 0; nu < Ndim_; ++nu){
      if(nu==mu) continue;

      mult_isigma(eta2, eta, mu, nu); //sigma_{mu,nu} mult

      U_nu.U = (*u_)[gf_->dir_slice(nu)];
      U_mu.U = (*u_)[gf_->dir_slice(mu)]; 

      //term 1 and 5 --------------------------------------
      Cmu_up.U  = stpl_->upper(*u_, mu,nu); //V_mu
      Cmu_up.U -= stpl_->lower(*u_, mu,nu); //V_+mu - V_-mu 
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = u_dag(Cmu_up,site) * v(eta2,s,site);
	  vect = u(U_mu,site) * vect;
	  vright.set(ff_->cslice(s,site), vect.getva());
	} 
      }
      external_prod(fce_tmp1.U, zeta, vright);
      fce_tmp2.U = fce_tmp1.U;

      //term 2 --------------------------------------
      Cnu_up.U = stpl_->upper(*u_, nu, mu);
      sf_up_[nu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(Cnu_up.U,site))* v(sf_up_[nu],s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      }
      sf_up_[nu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      }
      external_prod(fce_tmp1.U, vleft, vright);
      fce_tmp2.U += fce_tmp1.U;

      //term 4 and 8 --------------------------------------
      sf_up_[mu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      }    
      sf_up_[mu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(Cmu_up.U,site))* v(sf_up_[mu],s,site);
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      }    
      external_prod(fce_tmp1.U, vleft, vright);
      fce_tmp2.U += fce_tmp1.U;
   
      //term 3 --------------------------------------
      sf_up_[nu]->setf(eta2);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      }        
      sf_up_[mu]->setf(vright); //U_nu(x+mu)*eta2(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      } //vright = U_mu(x)*U_nu(x+mu)*eta2(x+mu+nu)  
     
      sf_up_[mu]->setf(zeta);
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      }        
      sf_up_[nu]->setf(vleft); //U_mu(x+nu)*zeta(x+mu+nu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_nu.U,site))* v(sf_up_[nu],s,site);
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      } //vright = U_nu(x)*U_mu(x+nu)*zeta(x+mu+nu)  
      external_prod(fce_tmp1.U, vleft, vright);
      fce_tmp2.U += fce_tmp1.U;

      //term 6 --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(Cnu_up.U,site))* v(eta2,s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      }       
      sf_dn_[nu]->setf(vright); // V_+nu (x-nu)*eta2(x-nu)
      shifted = sf_dn_[nu]->getva();

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(U_nu.U,site))* v(zeta,s,site);
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      }   
      sf_dn_[nu]->setf(vleft);
      external_prod(fce_tmp1.U, (Field)(sf_dn_[nu]->getva()), shifted);
      fce_tmp2.U -= fce_tmp1.U;   

      // term 7  --------------------------------------
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u_dag(U_nu.U,site))* v(eta2,s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      }   
      sf_dn_[nu]->setf(vright);
      shifted = sf_dn_[nu]->getva(); // shifted = Udag_nu(x-nu)*eta2(x-nu)
      sf_up_[mu]->setf(shifted);// Udag_nu(x+mu-nu)*eta2(x+mu-nu)

      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vright.set(ff_->cslice(s,site), vect.getva());
	}
      } //vright = U_mu(x)*Udag_nu(x+mu-nu)*eta2(x+mu-nu)
     
      sf_up_[mu]->setf(zeta); //zeta(x+mu)
      for(int site = 0; site<Nvol_; ++site){
	for(int s = 0; s < Nd; ++s) {
	  vect = (u(U_mu.U,site))* v(sf_up_[mu],s,site);
	  vect = (u_dag(U_nu.U,site))*vect;
	  vleft.set(ff_->cslice(s,site), vect.getva());
	}
      }       
      sf_dn_[nu]->setf(vleft);
      external_prod(fce_tmp1.U, (Field)(sf_dn_[nu]->getva()), vright);
      fce_tmp2.U -= fce_tmp1.U;   
     
      fce_tmp2.U *= - Dw->getKappa() * csw_ / 8.0; 

      for(int site = 0; site<Nvol_; ++site){
	force.U.add(force.Format.cslice(0,site,mu), 
		    anti_hermite(u(fce_tmp2.U,site)));

      }
    }
  }

  return force.U;
}

const vector<int> Dirac_Clover::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}
