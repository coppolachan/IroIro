#include "baryonCorrelator.hpp"
#include "Main/Geometry/siteIndex.hpp"
#include "Tools/gammaMatrices.hpp"
#include "Tools/colorEpsilon.hpp"
#include "Communicator/comm_io.hpp"
#include "include/field.h"
#include <complex>
using namespace std;
using namespace GammaMatrices;
using namespace ColorUtils;

complex<double> 
BaryonCorrelator::spTrTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			 const Gamma& G,int site,int i,int j)const{
  vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  epsilon_(fi,ei,i);
  epsilon_(fj,ej,j);

  double tr1r = 0.0, tr1i = 0.0;
  double tr2r = 0.0, tr2i = 0.0;

  for(int s=0; s<Nd_; ++s){
    tr1r += S1[Nc_*s+ei[0]][fmt_.index_r(ej[0],s,site)];
    tr1i += S1[Nc_*s+ei[0]][fmt_.index_i(ej[0],s,site)];
    
    GammaResult g = G(s);
    for(int sp=0; sp<Nd_; ++sp){
      GammaResult gp = G(sp);
      double facr = g.facr*gp.facr -g.faci*gp.faci;
      double faci = g.facr*gp.faci +g.faci*gp.facr;

      double tmpr 
	= S2[Nc_*g.spinor+ei[1]][fmt_.index_r(ej[1],gp.spinor,site)]*facr
	- S2[Nc_*g.spinor+ei[1]][fmt_.index_i(ej[1],gp.spinor,site)]*faci;

      double tmpi 
	= S2[Nc_*g+ei[1]][fmt_.index_r(ej[1],gp.spinor,site)]*faci
	+ S2[Nc_*g+ei[1]][fmt_.index_i(ej[1],gp.spinor,site)]*facr;

      tr2r += tmpr*S3[Nc_*sp+ei[2]][fmt_.index_r(ej[2],s)]
	     -tmpi*S3[Nc_*sp+ei[2]][fmt_.index_i(ej[2],s)];

      tr2i += tmpr*S3[Nc_*sp+ei[2]][fmt_.index_i(ej[2],s)]
             +tmpi*S3[Nc_*sp+ei[2]][fmt_.index_r(ej[2],s)];
    }
  }
  return complex<double>(tr1r*tr2r -tr1i*tr1i, tr1r*tr2i +tr1i*tr1r)*fi*fj; 
}

complex<double> 
BaryonCorrelator::spTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
		       const Gamma& G,int site,int i,int j)const{
  vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  epsilon_(fi,ei,i);
  epsilon_(fj,ej,j);

  double tr_r = 0.0, tr_i = 0.0;

  for(int sa=0; sa<Nd_; ++sa){
    for(int sb=0; sb<Nd_; ++sb){
      double S1r = S1[Nc_*sa+ei[0]][fmt_.index_r(ej[0],sb,site)];
      double S1i = S1[Nc_*sa+ei[0]][fmt_.index_i(ej[0],sb,site)];
      
      double S2S3r = 0.0, S2S3i = 0.0;

      GammaResult g = G(sb);
      for(int sc=0; sc<Nd_; ++sc){
	GammaResult gp = G(sc);
	double facr = g.facr*gp.facr -g.faci*gp.faci;
	double faci = g.facr*gp.faci +g.faci*gp.facr;

	double tmpr 
	  = S2[Nc_*g.spinor+ei[1]][fmt_.index_r(ej[1],gp.spinor,site)]*facr
	  - S2[Nc_*g.spinor+ei[1]][fmt_.index_i(ej[1],gp.spinor,site)]*faci;

	double tmpi 
	  = S2[Nc_*g+ei[1]][fmt_.index_r(ej[1],gp.spinor,site)]*faci
	  + S2[Nc_*g+ei[1]][fmt_.index_i(ej[1],gp.spinor,site)]*facr;

	S2S3r += tmpr*S3[Nc_*sc+ei[2]][fmt_.index_r(ej[2],sa)]
	        -tmpi*S3[Nc_*sc+ei[2]][fmt_.index_i(ej[2],sa)];

	S2S3i += tmpr*S3[Nc_*sc+ei[2]][fmt_.index_i(ej[2],sa)]
    	        +tmpi*S3[Nc_*sc+ei[2]][fmt_.index_r(ej[2],sa)];
      }
      tr_r += S1r*S2S3r -S1i*S2S3i;
      tr_i += S1r*S2S3i +S1i*S2S3r;
    }
  }
  return complex<double>(tr_r,tr_i)*fi*fj;
}

void global_correl(std::vector<complex<double> >& correl,
		   const std::vector<complex<double> >& correl_local){

  vector<complex<double> > correl_tmp(correl.size(),0.0);
  for(int t=0; t<Nt_; ++t) 
    correl_tmp[SiteIndex::instance()->global_t(t)] = correl_local[t];
  
  for(int t=0; t<Lt; ++t)
    correl[t] = complex<double>(comm_->reduce_sum(correl_tmp[t].real()),
				comm_->reduce_sum(correl_tmp[t].imag()));
}

const vector<complex<double> > octet(const prop_t& S1,const prop_t& S2){
  CGamma5 G;
  vector<complex<double> > correl_local(Nt_,0.0);
  
  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index

    for(int i=0; i<Nep_; ++i){
      for(int j=0; j<Nep_; ++j){
	correl_local[t] += spTrTr(S1,S1,S2,G,site,i,j);
	correl_local[t] += spTr(  S1,S1,S2,G,site,i,j);
      }
    }  
  }
  vector<complex<double> > correl_global(CommonPrms::instance()->Lt(),0.0);
  global_correl(correl_global,correl_local);
  return correl_global;
}

const vector<complex<double> > deculpet_del(const prop_t& S,int k){
  vector<complex<double> > correl_local(Nt_,0.0);
  Gamma* Gk; 
  if(k==XDIR) Gk = new CGamma1;
  if(k==YDIR) Gk = new CGamma2;
  if(k==ZDIR) Gk = new CGamma3;
  
  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index

    for(int i=0; i<Nep_; ++i){
      for(int j=0; j<Nep_; ++j){
	correl_local[t] += spTrTr(S1,S1,S2,Gk,site,i,j);
	correl_local[t] += spTr(  S1,S1,S2,Gk,site,i,j);
      }
    }  
  }
  delete Gk;

  vector<complex<double> > correl_global(CommonPrms::instance()->Lt(),0.0);
  global_correl(correl_global,correl_local);
  return correl_global;
}

const vector<complex<double> > deculpet_sgm(const prop_t& S1,const prop_t& S2,
					    int k){
  vector<complex<double> > correl_local(Nt_,0.0);
  Gamma* Gk; 
  if(k==XDIR) Gk = new CGamma1;
  if(k==YDIR) Gk = new CGamma2;
  if(k==ZDIR) Gk = new CGamma3;

  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);//get t from site index

    for(int i=0; i<Nep_; ++i){
      for(int j=0; j<Nep_; ++j){
	correl_local[t] += spTrTr(S1,S1,S2,*Gk,site,i,j);
	correl_local[t] += spTr(  S1,S1,S2,*Gk,site,i,j);
      }
    }  
  }
  delete Gk;

  vector<complex<double> > correl_global(CommonPrms::instance()->Lt(),0.0);
  global_correl(correl_global,correl_local);
  return correl_global;
}

const vector<complex<double> > nucleon(){ return octet(Sl_,Sl_);}
const vector<complex<double> > sigma8() { return octet(Sl_,Sh_);}
const vector<complex<double> > xi8()    { return octet(Sh_,Sl_);}
const vector<complex<double> > lambda() {}

const vector<complex<double> > delta(int k)  {return decouplet_del(Sl_,k);}
const vector<complex<double> > omega(int k)  {return decouplet_del(Sh_,k);}
const vector<complex<double> > sigma10(int k){return decouplet_sgm(Sl_,Sh_,k);}
const vector<complex<double> > xi10(int k)   {return decouplet_sgm(Sh_,Sl_,k);}

