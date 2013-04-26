#include "baryonCorrelator.hpp"
#include "Communicator/comm_io.hpp"
#include "include/field.h"
#include <complex>
using namespace std;
using namespace GammaMatrices;
using namespace ColorUtils;

complex<double> 
BaryonCorrelator::spTrTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			 const Gamma& G,int site,UpDn ud)const{
  int s_ini = ud*Nd_/2;
  int s_fin = (1+ud)*Nd_/2;

  vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  complex<double> correl(0.0,0.0);
  
  for(int i=0; i<Nep_; ++i){
    epsilon_.get_elm(fi,ei,i);

    for(int j=0; j<Nep_; ++j){
      epsilon_.get_elm(fj,ej,j);

      complex<double> eps(fi*fj,0.0);
      complex<double> Tr1(0.0,0.0);
      complex<double> Tr2(0.0,0.0);

      for(int sa=s_ini; sa<s_fin; ++sa){ /* contraction of a half-spinor*/
	
	Tr1 += complex<double>(S1[Nc_*sa+ei[0]][fmt_.index_r(ej[0],sa,site)],
			       S1[Nc_*sa+ei[0]][fmt_.index_i(ej[0],sa,site)]);
	GammaResult Ga = G(sa);

	for(int sb=0; sb<Nd_; ++sb){
	  GammaResult Gb = G(sb);
	  complex<double> fac(Ga.facr*Gb.facr -Ga.faci*Gb.faci,
			      Ga.facr*Gb.faci +Ga.faci*Gb.facr);
	  /*!@brief CgS2Cg should be transversed for spinor idx. */
	  complex<double> 
	    CgS2Cg(S2[Nc_*Gb.spn+ei[1]][fmt_.index_r(ej[1],Ga.spn,site)],
		   S2[Nc_*Gb.spn+ei[1]][fmt_.index_i(ej[1],Ga.spn,site)]);
	  CgS2Cg *= fac;
	  
	  Tr2 += CgS2Cg
	    *complex<double>(S3[Nc_*sb+ei[2]][fmt_.index_r(ej[2],sa,site)],
			     S3[Nc_*sb+ei[2]][fmt_.index_i(ej[2],sa,site)]);
	}
      }
      correl += eps*Tr1*Tr2;
    }
  }
  return correl;
}

complex<double> 
BaryonCorrelator::spTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
		       const Gamma& G,int site,UpDn ud)const{
  int s_ini = ud*Nd_/2;
  int s_fin = (1+ud)*Nd_/2;

  vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  complex<double> correl(0.0,0.0);
  
  for(int i=0; i<Nep_; ++i){
    epsilon_.get_elm(fi,ei,i);

    for(int j=0; j<Nep_; ++j){
      epsilon_.get_elm(fj,ej,j);

      complex<double> eps(fi*fj,0.0);
      complex<double> Tr(0.0,0.0);

      for(int sa=s_ini; sa<s_fin; ++sa){ /* contraction of a half-spinor*/
	for(int sb=0; sb<Nd_; ++sb){
	  complex<double> 
	    S1ab(S1[Nc_*sa+ei[0]][fmt_.index_r(ej[0],sb,site)],
		 S1[Nc_*sa+ei[0]][fmt_.index_i(ej[0],sb,site)]);
	  
	  complex<double> S3S2ba(0.0,0.0);
	  GammaResult Gb = G(sb);
	  
	  for(int sc=0; sc<Nd_; ++sc){
	    GammaResult Gc = G(sc);
	    complex<double> fac(Gb.facr*Gc.facr -Gb.faci*Gc.faci,
				Gb.facr*Gc.faci +Gb.faci*Gc.facr);
	    /*!@brief CgS3Cg should be transversed for spinor idx. */
	    complex<double> 
	      CgS3Cg(S3[Nc_*Gc.spn+ei[1]][fmt_.index_r(ej[1],Gb.spn,site)],
		     S3[Nc_*Gc.spn+ei[1]][fmt_.index_i(ej[1],Gb.spn,site)]);
	    CgS3Cg *= fac;
	    
	    S3S2ba += CgS3Cg
	      *complex<double>(S2[Nc_*sc+ei[2]][fmt_.index_r(ej[2],sa,site)],
			       S2[Nc_*sc+ei[2]][fmt_.index_i(ej[2],sa,site)]);
	  }
	  Tr += S1ab*S3S2ba;
	}
      }
      correl += eps*Tr;
    }
  }
  return correl;
}

correl_t BaryonCorrelator::
global_correl(const vector<complex<double> >& correl_local){

  int Lt = CommonPrms::instance()->Lt();
  correl_t correl(Lt,0.0);
  correl_t correl_tmp(Lt,0.0);

  for(int t=0; t<Nt_; ++t) 
    correl_tmp[SiteIndex::instance()->global_t(t)] = correl_local[t];
  
  for(int t=0; t<Lt; ++t)
    correl[t] = complex<double>(comm_->reduce_sum(correl_tmp[t].real()),
				comm_->reduce_sum(correl_tmp[t].imag()));
  return correl;  
}

/////// octet correlators //////
const correl_t 
BaryonCorrelator::octet(const prop_t& S1,const prop_t& S2,UpDn ud){
  CGamma5 G;

  vector<complex<double> > correl_local(Nt_,0.0);

  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);
    
    correl_local[t] += spTrTr(S1,S1,S2,G,site,ud);
    correl_local[t] += spTr(S1,S1,S2,G,site,ud);
  }
  return global_correl(correl_local);
}

const correl_t BaryonCorrelator::nucleon(UpDn ud){return octet(Sl_,Sl_,ud);}
const correl_t BaryonCorrelator::sigma8(UpDn ud){ return octet(Sl_,Sh_,ud);}
const correl_t BaryonCorrelator::xi8(UpDn ud){    return octet(Sh_,Sl_,ud);}

/////// singlet correlator //////
const correl_t BaryonCorrelator::lambda(UpDn ud){
  CGamma5 G;
  complex<double> r4(4.0),r2(2.0);
  vector<complex<double> > correl_local(Nt_,0.0);

  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);
    correl_local[t] += r2*spTrTr(Sl_,Sl_,Sh_,G,site,ud);
    correl_local[t] += r4*spTrTr(Sh_,Sl_,Sl_,G,site,ud);
    correl_local[t] += r4*spTr(  Sl_,Sh_,Sl_,G,site,ud);
    correl_local[t] += r4*spTr(  Sh_,Sl_,Sl_,G,site,ud);
    correl_local[t] -= r2*spTr(  Sl_,Sl_,Sh_,G,site,ud);
  }
  return global_correl(correl_local);
}

/////// decuplet correlators //////
const correl_t 
BaryonCorrelator::decuplet_del(const prop_t& S,site_dir k,UpDn ud){
  Gamma* Cg;
  if(k==XDIR) Cg = new CGamma1;
  if(k==YDIR) Cg = new CGamma2;
  if(k==ZDIR) Cg = new CGamma3;

  vector<complex<double> > correl_local(Nt_,0.0);
  complex<double> r2(2.0);

  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);

    correl_local[t] +=  spTrTr(S,S,S,*Cg,site,ud);
    correl_local[t] += r2*spTr(S,S,S,*Cg,site,ud);
  }
  delete Cg;
  return global_correl(correl_local);
}

const correl_t BaryonCorrelator::decuplet_sgm(const prop_t& S1,const prop_t& S2,
					      site_dir k,UpDn ud){
  Gamma* Cg;
  if(k==XDIR) Cg = new CGamma1;
  if(k==YDIR) Cg = new CGamma2;
  if(k==ZDIR) Cg = new CGamma3;

  vector<complex<double> > correl_local(Nt_,0.0);
  complex<double> r2(2.0);

  for(int site=0; site<Nvol_; ++site){
    int t = SiteIndex::instance()->c_t(site);

    correl_local[t] +=    spTrTr(S2,S1,S1,*Cg,site,ud);
    correl_local[t] += r2*spTrTr(S1,S1,S2,*Cg,site,ud);
    correl_local[t] += r2*spTr(  S1,S1,S2,*Cg,site,ud);
    correl_local[t] -= r2*spTr(  S1,S2,S1,*Cg,site,ud);
  }
  delete Cg;
  return global_correl(correl_local);
}

const correl_t BaryonCorrelator::delta(site_dir k,UpDn ud){
  return decuplet_del(Sl_,k,ud);
}
const correl_t BaryonCorrelator::omega(site_dir k,UpDn ud){
  return decuplet_del(Sh_,k,ud);
}
const correl_t BaryonCorrelator::sigma10(site_dir k,UpDn ud){
  return decuplet_sgm(Sl_,Sh_,k,ud);
}
const correl_t BaryonCorrelator::xi10(site_dir k,UpDn ud){
  return decuplet_sgm(Sh_,Sl_,k,ud);
}

