//----------------------------------------------------------------------
// dirac_wilson_Brillouin.cpp
//----------------------------------------------------------------------
#include "dirac_wilson_Brillouin.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Tools/randNum_MP.h"
#include "Communicator/comm_io.hpp"
#include "Fields/field_expressions.hpp"

using namespace SUNvecUtils;
using namespace std;

/////////////////////////////////////////////////////////////////////////
void Dirac_Wilson_Brillouin::get_RandGauss(valarray<double>& phi,
					   const RandNum& rng)const{
  MPrand::mp_get_gauss(phi,rng,SiteIndex::instance()->get_gsite(),ff_);
}

int Dirac_Wilson_Brillouin::sgm(int mu,int nu,int rho)const{
  for(int sigma=0; sigma<NDIM_; ++sigma){
    if(sigma !=mu && sigma !=nu && sigma !=rho) return sigma;
  }
}

const Field (Dirac_Wilson_Brillouin::*Dirac_Wilson_Brillouin::gm[])
(const Field&) const = {&Dirac_Wilson_Brillouin::gamma_x,
			&Dirac_Wilson_Brillouin::gamma_y,
			&Dirac_Wilson_Brillouin::gamma_z,
			&Dirac_Wilson_Brillouin::gamma_t,};

const Field Dirac_Wilson_Brillouin::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

const Field Dirac_Wilson_Brillouin::mult_H(const Field& f)const{ 
  return gamma5(mult(f));
}

/*!@brief U_{x,\mu}\psi_{x+\hat{\mu}} */

inline const Field Dirac_Wilson_Brillouin::sft_p(const Field& f,int dir) const{
  Field fp(fsize_);
  
  /// boundary part ///
  int is = 0;
  int Xb = 0;
  int Nbdry = slsize(dir); // ?
  double vbd[Nin_*Nbdry]; /*!< @brief data on the lower slice */   
  for(int k=0; k<Nbdry; ++k){
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xsl(Xb,k,dir)));
    // int index(int i,int site,int ex =0) <-- internal degree of freedom at site
    for(int s=0; s<ND_; ++s){
      for(int c=0; c<NC_; ++c){
	vbd[sr(s,c)+is] = v[sr(s,c)];
	vbd[si(s,c)+is] = v[si(s,c)];
      }
    }
    is += Nin_;
  }
  double vbc[Nin_*Nbdry]; /*!< @brief data on the upper slice */
  comm_->transfer_fw(vbc,vbd,Nin_*Nbdry,dir);

  is = 0;
  Xb = SiteIndex::instance()->Bdir(dir); // ?
  Nbdry = slsize(dir);

  for(int k=0; k<Nbdry; ++k){  /*!< @brief calc on the upper boundary */   

    int xc = xsl(Xb,k,dir); 
    const double* U = const_cast<Field*>(u_)->getaddr(gf_.index(0,xc,dir));
    double* res = fp.getaddr(ff_.index(0,xc));
    
    for(int s=0; s<ND_; ++s){
      for(int c=0; c<NC_; ++c){
	double vr = 0.0, vi = 0.0;
      
	for(int c1=0; c1<NC_; ++c1){
	  vr += U[re(c,c1)]*vbc[sr(s,c1)+is] -U[im(c,c1)]*vbc[si(s,c1)+is];
	  vi += U[im(c,c1)]*vbc[sr(s,c1)+is] +U[re(c,c1)]*vbc[si(s,c1)+is];
	}
	res[sr(s,c)] += vr;
	res[si(s,c)] += vi;
      }
    }
    is += Nin_;
  }
  /// bulk part ///
  for(int x=0; x<Xb; ++x){
    int Nslice = slsize(dir);
    for(int k=0; k<Nslice; ++k){   
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xsl(x+1,k,dir)));
      const double* U = const_cast<Field*>(u_)->getaddr(gf_.index(0,xsl(x,k,dir),dir));
      double* res = fp.getaddr(ff_.index(0,xsl(x,k,dir)));

      for(int s=0; s<ND_; ++s){
	for(int c=0; c<NC_; ++c){
	  double vr = 0.0, vi = 0.0;
	  
	  for(int c1=0; c1<NC_; ++c1){
	    vr += U[re(c,c1)]*v[sr(s,c1)] -U[im(c,c1)]*v[si(s,c1)];
	    vi += U[im(c,c1)]*v[sr(s,c1)] +U[re(c,c1)]*v[si(s,c1)];
	  }
	  res[sr(s,c)] += vr; 
	  res[si(s,c)] += vi;
	}
      }
    }
  }
  return fp;
}

/*! @brief  U_{x-\hat{\mu},\mu}^\dagger \psi_{x-\hat{\mu}} */

inline const Field Dirac_Wilson_Brillouin::sft_m(const Field& f,int dir)const{
  Field fm(fsize_);

  /// boundary part ///
  int is = 0;
  int Xb = SiteIndex::instance()->Bdir(dir);
  int Nbdry = slsize(dir);
  double vbd[Nin_*Nbdry]; /*!< @brief data on the upper boundary */
  for(int k=0; k<Nbdry; ++k){
    int xc = xsl(Xb,k,dir);
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xc));
    const double* U = const_cast<Field*>(u_)->getaddr(gf_.index(0,xc,dir));

    for(int s=0; s<ND_; ++s){    
      for(int c=0; c<NC_; ++c){
	vbd[sr(s,c)+is] = 0.0;  
	vbd[si(s,c)+is] = 0.0;
	
	for(int c1=0; c1<NC_; ++c1){
	  vbd[sr(s,c)+is] += U[re(c1,c)]*v[sr(s,c1)] +U[im(c1,c)]*v[si(s,c1)];
	  vbd[si(s,c)+is] -= U[im(c1,c)]*v[sr(s,c1)] -U[re(c1,c)]*v[si(s,c1)];
	}
      }
    }
    is += Nin_;
  }
  double vbc[Nin_*Nbdry];  //Copy vbd from backward processor
  comm_->transfer_bk(vbc,vbd,Nin_*Nbdry,dir);
  is = 0;
  Nbdry = slsize(dir);
  
  for(int k=0; k<Nbdry; ++k){
    double* res = fm.getaddr(ff_.index(0,xsl(0,k,dir)));
    
    for(int s=0; s<ND_; ++s){  
      for(int c=0; c<NC_; ++c){  
	res[sr(s,c)] += vbc[sr(s,c)+is];  
	res[si(s,c)] += vbc[si(s,c)+is];
      }
    }
    is += Nin_;
  }

  /// bulk part ///
  for(int x=1; x<Xb+1; ++x){
    int Nslice = slsize(dir);
    for(int k=0; k<Nslice; ++k){
      int xm = xsl(x-1,k,dir);
      const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xm));
      const double* U = const_cast<Field*>(u_)->getaddr(gf_.index(0,xm,dir));
      double* res = fm.getaddr(ff_.index(0,xsl(x,k,dir)));

      for(int s=0; s<ND_; ++s){
	for(int c=0; c<NC_; ++c){
	  double vr = 0.0, vi = 0.0; 
	
	  for(int c1=0; c1<NC_; ++c1){
	    vr += U[re(c1,c)]*v[sr(s,c1)] +U[im(c1,c)]*v[si(s,c1)];
	    vi -= U[im(c1,c)]*v[sr(s,c1)] -U[re(c1,c)]*v[si(s,c1)];
	  }
	  res[sr(s,c)] += vr;
	  res[si(s,c)] += vi;
	}
      }
    }
  }
  return fm;
}

inline const Field Dirac_Wilson_Brillouin::del(const Field& f,int dir,double a)const{
  using namespace FieldExpression;
  Field w = sft_p(f,dir);
  w -= sft_m(f,dir);
  w *= a;
  return w;
}

inline const Field Dirac_Wilson_Brillouin::lap(const Field& f,int dir,double a)const{
  using  namespace FieldExpression;
  Field w = sft_p(f,dir);
  w += sft_m(f,dir);
  w *= a;
  return w;
}

const Field Dirac_Wilson_Brillouin::mult(const Field& f)const{
  using namespace FieldExpression;
  Field w = mult_del(f);
  w -= 0.5*mult_lap(f);
  w *= kbr_;
  w += f;
  return w;
}

const Field Dirac_Wilson_Brillouin::mult_del(const Field& f)const{
  using namespace FieldExpression;
  Field v(fsize_);
  //===preconcitiong===================================
  Field f0 = lap(f,0,0.25);  Field f1 = lap(f,1,0.25);
  Field f2 = lap(f,2,0.25);  Field f3 = lap(f,3,0.25);
  f0 += 4.0*f; f1 += 4.0*f;  f2 += 4.0*f; f3 += 4.0*f;
  //===================================================
  for(int mu = 0; mu<NDIM_; ++mu){
    Field w(fsize_);
    Field xi3(fsize_),eta1(fsize_);
    //===preconcitiong==========================================
    Field f0d = del(f0,mu,1.0); Field f1d = del(f1,mu,1.0);
    Field f2d = del(f2,mu,1.0); Field f3d = del(f3,mu,1.0);
    eta1 = del(f,mu,1.0);
    Field f0l = lap(eta1,0,0.25); Field f1l = lap(eta1,1,0.25);
    Field f2l = lap(eta1,2,0.25); Field f3l = lap(eta1,3,0.25);
    //==========================================================
    for(int nu=0; nu<NDIM_; ++nu){
      Field xi2(fsize_),eta3(fsize_);
      if(nu !=mu){
	for(int rho=0; rho<NDIM_; ++rho){
	  Field xi1(fsize_),eta2(fsize_);
	  if(rho !=nu && rho != mu){
	    int dir = sgm(mu,nu,rho);
	    //===accerate improvement==============================
	    if(dir == 0)     {xi1 += f0; eta2 += f0d; eta2 += f0l;}
	    else if(dir == 1){xi1 += f1; eta2 += f1d; eta2 += f1l;}
	    else if(dir == 2){xi1 += f2; eta2 += f2d; eta2 += f2l;}
	    else if(dir == 3){xi1 += f3; eta2 += f3d; eta2 += f3l;}
	    //=====================================================
	    xi2 += lap(xi1,rho,1.0/3);
	    eta3 += lap(eta2,rho,1.0/3);
	  }//rho if
	}// rho loop
	xi2 += 16.0*f;
	eta3 += del(xi2,mu,1.0);
	xi3 += lap(xi2,nu,0.5);
	w += lap(eta3,nu,0.5/432);
      }// nu if
    }// nu loop
    xi3 += 64.0*f;
    w += del(xi3,mu,1.0/432);
    v += (this->*gm[mu])(w);
  } //mu loop
  return v;
}

const Field Dirac_Wilson_Brillouin::mult_lap(const Field& f)const{
  using namespace FieldExpression;
  Field w(fsize_);
  //===operations at a most inside loop====
  Field f0 = lap(f,0,0.25);
  Field f1 = lap(f,1,0.25);
  Field f2 = lap(f,2,0.25);
  Field f3 = lap(f,3,0.25);
  f0 += 2.0*f; f1 +=2.0*f;
  f2 += 2.0*f; f3 +=2.0*f;
  Field f01=lap(f0,1,1.0/3.0);Field f02=lap(f0,2,1.0/3.0);Field f03 = lap(f0,3,1.0/3.0);
  Field f10=lap(f1,0,1.0/3.0);Field f12=lap(f1,2,1.0/3.0);Field f13 = lap(f1,3,1.0/3.0);
  Field f20=lap(f2,0,1.0/3.0);Field f21=lap(f2,1,1.0/3.0);Field f23 = lap(f2,3,1.0/3.0);
  Field f30=lap(f3,0,1.0/3.0);Field f31=lap(f3,1,1.0/3.0);Field f32 = lap(f3,2,1.0/3.0);
  //=======================================
  for(int mu=0; mu<NDIM_; ++mu){
    Field fn(fsize_);
    for(int nu=0; nu<NDIM_; ++nu){
      Field fr(fsize_);
      if(nu !=mu){
	for(int rho=0; rho<NDIM_; ++rho){
	  //Field ft(fsize_);
	  if(rho !=nu && rho != mu){
	    int dir = sgm(mu,nu,rho);
	    //===accerate improvement====
	    if(dir == 0){
	      if(rho==1) fr += f01;
	      else if(rho==2) fr += f02; 
	      else if(rho==3) fr += f03;
	    }
	    else if(dir == 1){
	      if(rho==0) fr += f10;
	      else if(rho==2) fr += f12;
	      else if(rho==3) fr += f13;
	    }
	    else if(dir == 2){
	      if(rho==0) fr += f20;
	      else if(rho==1) fr += f21;
	      else if(rho==3) fr += f23;
	    }
	    else if(dir == 3){
	      if(rho==0) fr += f30;
	      else if(rho==1) fr += f31;
	      else if(rho==2) fr += f32;
	    }
	    //===========================
	    //fr += lap(ft,rho,1.0/3);
	  }// rho if
	}// rho looop
	fr += 4.0*f;
	fn += lap(fr,nu,0.5); 
      }//nu if
    }// nu loop
    fn += 8.0*f;
    w += lap(fn,mu,1.0/64);
  }// mu loop
  return w;
}

// Dirac representation 
const Field Dirac_Wilson_Brillouin::gamma5(const Field& f) const{
  Field w(fsize_);

  for(int site=0; site<Nvol_*Nin_; site+=Nin_){
    const double* ft = const_cast<Field&>(f).getaddr(site);
    double* res = w.getaddr(site);    

    for(int c=0; c<NC_; ++c){
      short int cc = 2*c;
      short int cc1 = 2*c+1;

      res[      cc] = ft[4*NC_+cc];  res[      cc1] = ft[4*NC_+cc1];
      res[2*NC_+cc] = ft[6*NC_+cc];  res[2*NC_+cc1] = ft[6*NC_+cc1];
      res[4*NC_+cc] = ft[      cc];  res[4*NC_+cc1] = ft[      cc1];
      res[6*NC_+cc] = ft[2*NC_+cc];  res[6*NC_+cc1] = ft[2*NC_+cc1];
    }
  }
  return w;
}

const Field Dirac_Wilson_Brillouin::gamma_x(const Field& f) const{
  Field w(fsize_);

  for(int site=0; site<Nvol_*Nin_; site+=Nin_){
    const double* ft = const_cast<Field&>(f).getaddr(site);
    double* res = w.getaddr(site);    

    for(int c=0; c<NC_; ++c){
      short int cc = 2*c;
      short int cc1 = 2*c+1;
      
      res[      cc] = ft[6*NC_+cc1];  res[      cc1] =-ft[6*NC_+cc];
      res[2*NC_+cc] = ft[4*NC_+cc1];  res[2*NC_+cc1] =-ft[4*NC_+cc];
      res[4*NC_+cc] =-ft[2*NC_+cc1];  res[4*NC_+cc1] = ft[2*NC_+cc];
      res[6*NC_+cc] =-ft[      cc1];  res[6*NC_+cc1] = ft[      cc];
    }
  }
  return w;
}

const Field Dirac_Wilson_Brillouin::gamma_y(const Field& f) const{
  Field w(fsize_);
  
  for(int site=0; site<Nvol_*Nin_; site+=Nin_){
    const double* ft = const_cast<Field&>(f).getaddr(site);
    double* res = w.getaddr(site);    

    for(int c=0; c<NC_; ++c){
      short int cc = 2*c;
      short int cc1 = 2*c+1;
      
      res[      cc] =-ft[6*NC_+cc];  res[      cc1] =-ft[6*NC_+cc1];
      res[2*NC_+cc] = ft[4*NC_+cc];  res[2*NC_+cc1] = ft[4*NC_+cc1];
      res[4*NC_+cc] = ft[2*NC_+cc];  res[4*NC_+cc1] = ft[2*NC_+cc1];
      res[6*NC_+cc] =-ft[      cc];  res[6*NC_+cc1] =-ft[      cc1];
    }
  }
  return w;
}

const Field Dirac_Wilson_Brillouin::gamma_z(const Field& f) const{
  Field w(fsize_);

  for(int site=0; site<Nvol_*Nin_; site+=Nin_){
    const double* ft = const_cast<Field&>(f).getaddr(site);
    double* res = w.getaddr(site);    

    for(int c=0; c <NC_; ++c){
      short int cc = 2*c;
      short int cc1 = 2*c+1;
      
      res[      cc] = ft[4*NC_+cc1];  res[      cc1] =-ft[4*NC_+cc];
      res[2*NC_+cc] =-ft[6*NC_+cc1];  res[2*NC_+cc1] = ft[6*NC_+cc];
      res[4*NC_+cc] =-ft[      cc1];  res[4*NC_+cc1] = ft[      cc];
      res[6*NC_+cc] = ft[2*NC_+cc1];  res[6*NC_+cc1] =-ft[2*NC_+cc];
    }
  }
  return w;
}

const Field Dirac_Wilson_Brillouin::gamma_t(const Field& f) const{
  Field w(f);
  for(int site=0; site<Nvol_; ++site){
    w.set(ff_.cslice(2,site), -f[ff_.cslice(2,site)]);
    w.set(ff_.cslice(3,site), -f[ff_.cslice(3,site)]);
  }
  return w;
}
