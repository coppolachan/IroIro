//----------------------------------------------------------------------
// dirac_wilson_Brillouin_Imp.cpp
//----------------------------------------------------------------------
#include "dirac_wilson_Brillouin_Imp.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Tools/randNum_MP.h"
#include "Fields/field_expressions.hpp"

using namespace SUNvecUtils;
using namespace std;

void Dirac_Wilson_Brillouin_Imp::get_RandGauss(valarray<double>& phi,
					       const RandNum& rng)const{
  MPrand::mp_get(phi,rng,SiteIndex::instance()->get_gsite(),ff_);
}

int Dirac_Wilson_Brillouin_Imp::sgm(int mu,int nu, int rho)const{
  for(int sigma=0; sigma<NDIM_; ++sigma){
    if(sigma !=mu && sigma !=nu && sigma !=rho) return sigma;
  }
}

const Field (Dirac_Wilson_Brillouin_Imp::*Dirac_Wilson_Brillouin_Imp::gm[])
(const Field&) const = {&Dirac_Wilson_Brillouin_Imp::gamma_x,
			&Dirac_Wilson_Brillouin_Imp::gamma_y,
			&Dirac_Wilson_Brillouin_Imp::gamma_z,
			&Dirac_Wilson_Brillouin_Imp::gamma_t,};

const Field Dirac_Wilson_Brillouin_Imp::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

/*!@brief U_{x,\mu}\psi_{x+\hat{\mu}} */

const Field Dirac_Wilson_Brillouin_Imp::sft_p(const Field& f,int dir) const{
  Field fp(fsize_);

  /// boundary part ///
  int is = 0;
  int Xb = 0;
  int Nbdry = slsize(dir);
  double vbd[Nin_*Nbdry]; /*!< @brief data on the lower slice */   
  for(int k=0; k<Nbdry; ++k){
    const double* v = const_cast<Field&>(f).getaddr(ff_.index(0,xsl(Xb,k,dir)));
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
  Xb = SiteIndex::instance()->Bdir(dir);
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

const Field Dirac_Wilson_Brillouin_Imp::sft_m(const Field& f,int dir)const{
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

const Field Dirac_Wilson_Brillouin_Imp::delg(const Field& f,int dir,double a,double b)const{
  using namespace FieldExpression;
  Field w = sft_p(f,dir);
  w -= sft_m(f,dir);
  w *= a;
  w += b*f;
  return (this->*gm[dir])(w);
}  

const Field Dirac_Wilson_Brillouin_Imp::lap(const Field& f,int dir,double a,double b)const{
  using  namespace FieldExpression;
  Field w = sft_p(f,dir);
  w += sft_m(f,dir);
  w *= a;
  w += b*f;
  return w;
}  

const Field Dirac_Wilson_Brillouin_Imp::lap_hop(const Field& f)const{
  using namespace FieldExpression;
  Field w(fsize_);

  for(int mu=0; mu<NDIM_; ++mu){
    Field fn(fsize_);
    for(int nu=0; nu<NDIM_; ++nu){
      if(nu !=mu){
	Field fr(fsize_);
	for(int rho=0; rho<NDIM_; ++rho){
	  if(rho !=nu && rho != mu){
	    Field ft(fsize_);
	    int dir = sgm(mu,nu,rho);

	    ft += lap(f,dir,0.25,2.0);

	    fr += lap(ft,rho,1.0/3,0.0);
	  }
	}
	fr += 4.0*f;
	fn += lap(fr,nu,0.5,0.0);
      }
    }
    fn += 8.0*f;
    w += lap(fn,mu,1.0/64,0.0);
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_Imp::del_iso(const Field& f)const{
  using namespace FieldExpression;
  Field w(fsize_);

  for(int mu=0; mu<NDIM_; ++mu){
    Field dn(fsize_);
    for(int nu=0; nu<NDIM_; ++nu){
      if(nu !=mu){
	Field dr(fsize_);
	for(int rho=0; rho<NDIM_; ++rho){
	  if(rho !=nu && rho != mu){
	    Field dt(fsize_);
	    int dir = sgm(mu,nu,rho);

	    dt += lap(f,dir,1.0/3,4.0);

	    dr += lap(dt,rho,1.0/2,0.0);
	  }
	}
	dr += 16.0*f;
	dn += lap(dr,nu,1.0,0.0);
      }
    }
    dn += 64.0*f;
    w += delg(dn,mu,1.0/432,0.0);
  }
  return w;
}

const Field Dirac_Wilson_Brillouin_Imp::mult(const Field& f)const{
  using namespace FieldExpression;

  Field fn = lap_hop(lap_hop(f));
  fn -= 7.50*lap_hop(f);
  fn *= 0.125;

  Field dn = del_iso(f);
  dn -= 1.0/6.0*(lap_hop(dn)-15.0/4.0*dn);

  dn += fn;
  dn *= kimp_;
  dn += f;
  return dn;
}

const Field Dirac_Wilson_Brillouin_Imp::gamma5(const Field& f) const{
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

const Field Dirac_Wilson_Brillouin_Imp::gamma_x(const Field& f) const{
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

const Field Dirac_Wilson_Brillouin_Imp::gamma_y(const Field& f) const{
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

const Field Dirac_Wilson_Brillouin_Imp::gamma_z(const Field& f) const{
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

const Field Dirac_Wilson_Brillouin_Imp::gamma_t(const Field& f) const{
  Field w(f);
  for(int site=0; site<Nvol_; ++site){
    w.set(ff_.cslice(2,site), -f[ff_.cslice(2,site)]);
    w.set(ff_.cslice(3,site), -f[ff_.cslice(3,site)]);
  }
  return w;
}
