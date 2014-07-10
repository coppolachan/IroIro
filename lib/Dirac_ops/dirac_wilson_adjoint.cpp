/*! @file dirac_wilson_adjoint.cpp
 * @brief implementation of the Dirac_Wilson_Adjoint class 
 Time-stamp: <2014-06-05 16:35:06 noaki>
 */
#include "dirac_wilson_adjoint.hpp"
#include "Tools/sunMatUtils.hpp"

#include "Tools/fieldUtils.hpp"
#include "Tools/sunVec.hpp"

#include <omp.h>

using namespace std;

/* first half of the mult */
void Dirac_Wilson_Adjoint::mult_p(Field& af,const Field& f,int dir) const{
  Field Dpf(ff_.size()),ft(ff_.size());    
  for(int c=0; c<NC_; ++c){
    lmd_ad2fd(ft,f,c);
    Dpf = 0.0;
    (Dw_.*(Dw_.mult_p[dir]))(Dpf,ft);
    lmd_fd2Uad(af,Dpf,dir,c);
  }
}

void Dirac_Wilson_Adjoint::mult_offdiag(Field& Df,const Field& f) const{
  Field tmp(ff_.size()),ft(ff_.size());

  for(int c=0; c<NC_; ++c){
    lmd_ad2fd(ft,f,c);
    for(int mu=0; mu<NDIM_; ++mu){
      tmp = 0.0;
      (Dw_.*(Dw_.mult_p[mu]))(tmp,ft);
      lmd_fd2Uad(Df,tmp,mu,c);
    }
    for(int mu=0; mu<NDIM_; ++mu){
      lmd_ad2Ufd(ft,f,mu,c);
      tmp = 0.0; 
      (Dw_.*(Dw_.mult_m[mu]))(tmp,ft);
      lmd_fd2ad(Df,tmp,c);
    }
  }
  Df *= -0.5*kpp_;
}

void Dirac_Wilson_Adjoint::mult_full(Field& Df, const Field& f) const{
  mult_offdiag(Df,f);
  Df += f; 
}

const Field Dirac_Wilson_Adjoint::mult(const Field& f) const{
  Field Df(af_.size());
  (this->*mult_core)(Df,f);
  return Df;
}

////////////////////////////////////////////////////////////////////
/*** conversion: adjoint -> fundamental ***/

double (Dirac_Wilson_Adjoint::*Dirac_Wilson_Adjoint::laf[])(double* fp)const 
= {&Dirac_Wilson_Adjoint::laf11r, &Dirac_Wilson_Adjoint::laf11i,
   &Dirac_Wilson_Adjoint::laf12r, &Dirac_Wilson_Adjoint::laf12i,
   &Dirac_Wilson_Adjoint::laf13r, &Dirac_Wilson_Adjoint::laf13i,
   &Dirac_Wilson_Adjoint::laf21r, &Dirac_Wilson_Adjoint::laf21i,
   &Dirac_Wilson_Adjoint::laf22r, &Dirac_Wilson_Adjoint::laf22i,
   &Dirac_Wilson_Adjoint::laf23r, &Dirac_Wilson_Adjoint::laf23i,
   &Dirac_Wilson_Adjoint::laf31r, &Dirac_Wilson_Adjoint::laf31i,
   &Dirac_Wilson_Adjoint::laf32r, &Dirac_Wilson_Adjoint::laf32i,
   &Dirac_Wilson_Adjoint::laf33r, &Dirac_Wilson_Adjoint::laf33i,};

/* (\lambda^a\psi^a)_{c,c0} : a fundamental field */
void Dirac_Wilson_Adjoint::lmd_ad2fd(Field& lf,const Field& f,int c0)const{

#pragma omp parallel
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;
    /*
    for(int site=is; site<is+ns; ++site){
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(af_.index_r(0,s,site));      
	for(int c=0; c<NC_; ++c){
	  lf.set(ff_.index_r(c,s,site),(this->*laf[re(c,c0)])(fp));
	  lf.set(ff_.index_i(c,s,site),(this->*laf[im(c,c0)])(fp));
	}
      }
    }
    */
    for(int site=is; site<is+ns; ++site){
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(af_.index_r(0,s,site));      
	double* lfp = lf.getaddr(ff_.index_r(0,s,site));
	for(int c=0; c<NC_; ++c){
	  lfp[r(c)] = (this->*laf[re(c,c0)])(fp);
	  lfp[i(c)] = (this->*laf[im(c,c0)])(fp);
	}
      }
    }
  }
}

/* (\lambda^a\psi^a U_mu)_{c,c0} : a fundamental field */
void Dirac_Wilson_Adjoint::lmd_ad2Ufd(Field& lfu,const Field& f,int mu,int c0)const{
#pragma omp parallel
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;
    /*
    for(int site=is; site<is+ns; ++site){
      double* up = const_cast<Field*>(u_)->getaddr(gf_.index(0,gm_[site],mu));
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(af_.index_r(0,s,site));
	for(int c=0; c<NC_; ++c){
	  
	  complex<double> x(0.0,0.0);
	  for(int c1=0; c1<NC_; ++c1){
	    complex<double> fcc1((this->*laf[re(c,c1)])(fp),
				 (this->*laf[im(c,c1)])(fp));
	    complex<double> uc1c0(up[re(c1,c0)],up[im(c1,c0)]);
	    x += fcc1*uc1c0;
	  }
	  lfu.set(ff_.index_r(c,s,site),x.real());
	  lfu.set(ff_.index_i(c,s,site),x.imag());
	}
      }
    }
    */
    for(int site=is; site<is+ns; ++site){
      double* up = const_cast<Field*>(u_)->getaddr(gf_.index(0,gm_[site],mu));
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(af_.index_r(0,s,site));
	double* lfp = lfu.getaddr(ff_.index_r(0,s,site)) ;

	for(int c=0; c<NC_; ++c){
	  lfp[r(c)] = 0.0;
	  lfp[i(c)] = 0.0;
	  for(int c1=0; c1<NC_; ++c1){
	    lfp[r(c)] += (this->*laf[re(c,c1)])(fp)*up[re(c1,c0)] -(this->*laf[im(c,c1)])(fp)*up[im(c1,c0)];
	    lfp[i(c)] += (this->*laf[im(c,c1)])(fp)*up[re(c1,c0)] +(this->*laf[re(c,c1)])(fp)*up[im(c1,c0)];
	  }
	}
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
/*** conversion: fundamental -> adjoint ***/
/* (\lambda^a\Psi)_{c0c0} */
void Dirac_Wilson_Adjoint::lfa(valarray<double>& h,double* fp,int c0)const {
  h = 0.0;
  if(c0==0){
    h[0] = fp[r(1)]; h[1] = fp[i(1)];
    h[2] = fp[i(1)]; h[3] =-fp[r(1)];
    h[4] = fp[r(0)]; h[5] = fp[i(0)];
    h[6] = fp[r(2)]; h[7] = fp[i(2)];
    h[8] = fp[i(2)]; h[9] =-fp[r(2)];
    h[14]= fp[r(0)]/sqrt(3.0); h[15]= fp[i(0)]/sqrt(3.0);
  }else if(c0==1){
    h[0] = fp[r(0)]; h[1] = fp[i(0)];
    h[2] =-fp[i(0)]; h[3] = fp[r(0)];
    h[4] =-fp[r(1)]; h[5] =-fp[i(1)];
    h[10]= fp[r(2)]; h[11]= fp[i(2)];
    h[12]= fp[i(2)]; h[13]=-fp[r(2)];
    h[14]= fp[r(1)]/sqrt(3.0); h[15]= fp[i(1)]/sqrt(3.0);
  }else if(c0==2){
    h[6] = fp[r(0)];  h[7] = fp[i(0)];
    h[8] =-fp[i(0)];  h[9] = fp[r(0)];
    h[10]= fp[r(1)];  h[11]= fp[i(1)];
    h[12]=-fp[i(1)];  h[13]= fp[r(1)];
    h[14]= -2.0/sqrt(3.0)*fp[r(2)];  
    h[15]= -2.0/sqrt(3.0)*fp[i(2)];
  }else{
    abort();
  }
}

/* (Udag\lambda^a\Psi)_c0c0 */
void Dirac_Wilson_Adjoint::lfUa(valarray<double>& h,double* fp,
				double* up,int c0)const {

  complex<double> u0(up[re(0,c0)],-up[im(0,c0)]);
  complex<double> u1(up[re(1,c0)],-up[im(1,c0)]);
  complex<double> u2(up[re(2,c0)],-up[im(2,c0)]);
  
  complex<double> f0(fp[r(0)],fp[i(0)]);
  complex<double> f1(fp[r(1)],fp[i(1)]);
  complex<double> f2(fp[r(2)],fp[i(2)]);

  complex<double> t0 = f1*u0;
  complex<double> t1 = f0*u1;

  h[0] = t0.real() +t1.real();   h[1] = t0.imag() +t1.imag();
  h[2] = t0.imag() -t1.imag();   h[3] =-t0.real() +t1.real();
  
  t0 = f2*u0;   t1 = f0*u2;

  h[6] = t0.real() +t1.real();   h[7] = t0.imag() +t1.imag();
  h[8] = t0.imag() -t1.imag();   h[9] =-t0.real() +t1.real();

  t0 = f2*u1;  t1 = f1*u2;

  h[10] = t0.real() +t1.real();  h[11] = t0.imag() +t1.imag();
  h[12] = t0.imag() -t1.imag();  h[13] =-t0.real() +t1.real();
  
  t0 = f0*u0;  t1 = f1*u1;

  h[ 4] = t0.real() -t1.real();  h[ 5] = t0.imag() -t1.imag();
  h[14] = t0.real() +t1.real();  h[15] = t0.imag() +t1.imag();
  
  t0 = f2*u2;
  h[14] -= 2.0*t0.real();  h[15] -= 2.0*t0.imag();
  h[14] /= sqrt(3.0);      h[15] /= sqrt(3.0);
}

/* (\lambda^a\Psi)_{c0c0}*/
void Dirac_Wilson_Adjoint::lmd_fd2ad(Field& af,const Field& f,int c0)const{

#pragma omp parallel 
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;

    valarray<double> h(af_.Nin());    
    for(int site=is; site<is+ns; ++site){
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(ff_.index_r(0,s,site));
	lfa(h,fp,c0);
	af.add(af_.cslice(s,site),h);
      }
    }
  }
}

/* (Udmu\lambda^a\Psi)_{c0c0}*/
void Dirac_Wilson_Adjoint::lmd_fd2Uad(Field& af,const Field& f,int mu,int c0)const{

#pragma omp parallel 
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;

    valarray<double> h(af_.Nin());
    for(int site=is; site<is+ns; ++site){
      double* up = const_cast<Field*>(u_)->getaddr(gf_.index(0,gp_[site],mu));
      for(int s=0; s<ND_; ++s){
	double* fp = const_cast<Field&>(f).getaddr(ff_.index_r(0,s,site));
	lfUa(h,fp,up,c0);
	af.add(af_.cslice(s,site),h);
      }
    }
  }
}

const Field Dirac_Wilson_Adjoint::gamma5(const Field& f)const{ 

  Field w(af_.size());
#pragma omp parallel 
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;

    for(int site=is; site<is+ns; ++site){
      dm_.gamma5core(w.getaddr(af_.index(0,site)),
		     const_cast<Field&>(f).getaddr(af_.index(0,site)));
    }
  }
  return w;
}

const Field Dirac_Wilson_Adjoint::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

void Dirac_Wilson_Adjoint::mkfrc(Field& fce,const Field& eta,
				 const Field& zeta,int mu)const{
  using namespace SUNmatUtils;

  Field xie(af_.size());
  mult_p(xie,eta,mu);

#pragma omp parallel 
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;
    //SUNmat f;
    valarray<double> f(gf_.Nin());
    for(int site=is; site<is+ns; ++site){
      f = 0.0;
      double* ztp = const_cast<Field&>(zeta).getaddr(af_.index(0,site));
      double* xep = xie.getaddr(af_.index(0,site));

      for(int a=0; a<NADJ_; ++a){
	for(int b=0; b<NADJ_; ++b){

	  if(su3str[a*NADJ_+b] != 0){
	    double prod =0.0;
	    for(int s=0; s<ND_; ++s)
	      prod += ztp[ar(b,s)]*xep[ar(a,s)] +ztp[ai(b,s)]*xep[ai(a,s)];
	    f += lmd_commutator_va(a,b)*prod;
	  }
	}
      }
      f *= 0.25;
      //fce.add(gf_.islice(gp_[site],mu),f.getva());
      fce.add(gf_.islice(gp_[site],mu),f);
    } 
  }
}

void Dirac_Wilson_Adjoint::
md_force_p(Field& fce,const Field& eta,const Field& zeta)const{
  for(int mu=0; mu<NDIM_; ++mu) mkfrc(fce,eta,zeta,mu);
}  

void Dirac_Wilson_Adjoint::
md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  Field zt5 = gamma5(zeta);
  Field et5 = gamma5(eta);
  for(int mu=0; mu<NDIM_; ++mu) mkfrc(fce,zt5,et5,mu);
}  

const Field Dirac_Wilson_Adjoint::md_force(const Field& eta,
					   const Field& zeta)const{
  Field fce(gf_.size());
  md_force_p(fce,eta,zeta);
  md_force_m(fce,eta,zeta);
  fce *= -kpp_;
  return fce;
}

/*
//////////////////////////////////////////////////////////////////
/// for test purpose
const Field Dirac_Wilson_Adjoint::fadU(const Field& f,int mu)const{
  using namespace FieldUtils;
  using namespace SUNvecUtils;

  Field af(af_.size());
  valarray<double> fa(af_.Nin());

  for(int c=0; c<NC_; ++c){
    FermionField fd(lmd_ad2Ufd(f,mu,c));
    Field ff(ff_.size());

    for(int site=0; site<Nvol_; ++site){
      for(int sp=0; sp<ND_; ++sp){
	ff.set(ff_.cslice(sp,site),
	       (mat_dag(GaugeField(*u_),site,mu)*vec(fd,sp,site)).getva());
      }
    }
    lmd_fd2ad(af,ff,c);
  }
  af *= 0.5;
  return af;
}

const Field Dirac_Wilson_Adjoint::Ufad(const Field& f,int mu)const{
  using namespace FieldUtils;
  using namespace SUNvecUtils;

  Field af(af_.size());
  valarray<double> fa(af_.Nin());
  
  for(int c=0; c<NC_; ++c){
    FermionField fd(lmd_ad2fd(f,c));
    Field ff(ff_.size());
    
    for(int site=0; site<Nvol_; ++site){    
      for(int sp=0; sp<ND_; ++sp){
	ff.set(ff_.cslice(sp,site),
	       (mat(GaugeField(*u_),site,mu)*vec(fd,sp,site)).getva());
      }	
    }
    lmd_fd2Uad(af,ff,mu,c);
  }
  af *= 0.5;
  return af;
}

void Dirac_Wilson_Adjoint::test_unitarity(const Field& f,int mu)const{
  Field uf_u = fadU(Ufad(f,mu),mu);
  Field u_fu = Ufad(fadU(f,mu),mu);  

  for(int site=0; site<Nvol_; ++site){
    for(int s=0; s<ND_; ++s){
      for(int c=0; c<NADJ_; ++c){
	CCIO::cout<<site<<" "<<s<<" "<<c<<" real: "
		  <<f[   af_.index_r(c,s,site)]<<" "
		  <<uf_u[af_.index_r(c,s,site)]<<" "
		  <<u_fu[af_.index_r(c,s,site)]<<"\n";
	CCIO::cout<<site<<" "<<s<<" "<<c<<" imag: "
		  <<f[   af_.index_i(c,s,site)]<<" "
		  <<uf_u[af_.index_i(c,s,site)]<<" "
		  <<u_fu[af_.index_i(c,s,site)]<<"\n";
      }
    }
  }
}
///
*/

