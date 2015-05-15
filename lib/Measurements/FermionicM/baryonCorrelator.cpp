#include "baryonCorrelator.hpp"
#include "Communicator/comm_io.hpp"
#include "include/field.h"
#include <complex>
using namespace std;
using namespace GammaMatrices;
using namespace ColorUtils;

complex<double> 
BaryonCorrelator::spTrTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			 const Gamma& G, int site,UpDn ud)const{
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
BaryonCorrelator::spTrTrOpt(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			    const Gamma& G,
			    complex<double> &r1,complex<double> &r2,
			    int site,UpDn ud)const{
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

      complex<double> Tr3(0.0,0.0);
      complex<double> Tr4(0.0,0.0);

      for(int sa=s_ini; sa<s_fin; ++sa){ /* contraction of a half-spinor*/
	
	Tr1 += complex<double>(S1[Nc_*sa+ei[0]][fmt_.index_r(ej[0],sa,site)],
			       S1[Nc_*sa+ei[0]][fmt_.index_i(ej[0],sa,site)]);

	Tr3 += complex<double>(S3[Nc_*sa+ei[0]][fmt_.index_r(ej[0],sa,site)],
			       S3[Nc_*sa+ei[0]][fmt_.index_i(ej[0],sa,site)]);

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
	  Tr4 += CgS2Cg
	    *complex<double>(S1[Nc_*sb+ei[2]][fmt_.index_r(ej[2],sa,site)],
			     S1[Nc_*sb+ei[2]][fmt_.index_i(ej[2],sa,site)]);

	}
      }
      correl += eps*r1*Tr1*Tr2;
      correl += eps*r2*Tr3*Tr4;

    }
  }
  return correl;
}


complex<double> 
BaryonCorrelator::spTr(const prop_t& S1,const prop_t& S2,const prop_t& S3,
		       const Gamma& G,int site,UpDn ud)const{
  int s_ini = ud*Nd_/2;
  int s_fin = (1+ud)*Nd_/2;

  //vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  complex<double> correl(0.0,0.0);
  
  for(int i=0; i<Nep_; ++i){
    //epsilon_.get_elm(fi,ei,i);
    fi = 1.0-2.0*(i/Nc_);

    for(int j=0; j<Nep_; ++j){
      //epsilon_.get_elm(fj,ej,j);
      fj = 1.0-2.0*(j/Nc_);


      //complex<double> eps(fi*fj,0.0);
      complex<double> eps(fi*fj,0.0);
      complex<double> Tr(0.0,0.0);


      int base_idx = fmt_.index_r(eps_[j*3],0,site);

      for(int sa=s_ini; sa<s_fin; ++sa){ /* contraction of a half-spinor*/
	int sa_idx = Nc_*sa+eps_[i*3];

	for(int sb=0; sb<Nd_; ++sb){
	  /*
	  complex<double> 
	    S1ab(S1[sa_idx][fmt_.index_r(ej[0],sb,site)],
		 S1[sa_idx][fmt_.index_i(ej[0],sb,site)]);
	  */
	  complex<double> 
	    S1ab(S1[sa_idx][base_idx+ 2*Nc_*sb],
		 S1[sa_idx][base_idx+ 2*Nc_*sb+1]);


	  complex<double> S3S2ba(0.0,0.0);
	  GammaResult Gb = G(sb);
	  
	  for(int sc=0; sc<Nd_; ++sc){
	    GammaResult Gc = G(sc);
	    //common factor that we do not want to recalculate for every site...
	    //fac is an Nd*Nd matrix
	    //instead of passing G pass the precalculated matrix
	    complex<double> fac(Gb.facr*Gc.facr -Gb.faci*Gc.faci,
	    			Gb.facr*Gc.faci +Gb.faci*Gc.facr);

	   
	    /*!@brief CgS3Cg should be transversed for spinor idx. */
	    complex<double> 
	      CgS3Cg(S3[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_r(eps_[j*3+1],Gb.spn,site)],
		     S3[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_i(eps_[j*3+1],Gb.spn,site)]);
	    CgS3Cg *= fac;
	    
	    S3S2ba += CgS3Cg
	      *complex<double>(S2[Nc_*sc+eps_[i*3+2]][fmt_.index_r(eps_[j*3+2],sa,site)],
			       S2[Nc_*sc+eps_[i*3+2]][fmt_.index_i(eps_[j*3+2],sa,site)]);
	  }
	  Tr += S1ab*S3S2ba;
	}
      }
      correl += eps*Tr;
    }
  }
  return correl;
}


complex<double> 
BaryonCorrelator::spTrOpt(const prop_t& S1,const prop_t& S2,const prop_t& S3,
			  const Gamma& G, 
			  complex<double> &r1,complex<double> &r2,
			  vector< complex<double > > &fac_matrix, 
			  int site,UpDn ud)const{
  int s_ini = ud*Nd_/2;
  int s_fin = (1+ud)*Nd_/2;

  //vector<int> ei(Nc_),ej(Nc_);
  double fi,fj; 
  complex<double> correl(0.0,0.0);
  
  for(int i=0; i<Nep_; ++i){
    fi = 1.0-2.0*(i/Nc_);

    for(int j=0; j<Nep_; ++j){
      fj = 1.0-2.0*(j/Nc_);

      complex<double> eps(fi*fj,0.0);
      complex<double> Tr1(0.0,0.0);
      complex<double> Tr2(0.0,0.0);


      int base_idx = fmt_.index_r(eps_[j*3],0,site);

      for(int sa=s_ini; sa<s_fin; ++sa){ /* contraction of a half-spinor*/
	int sa_idx = Nc_*sa+eps_[i*3];

	for(int sb=0; sb<Nd_; ++sb){
	  /*
	  complex<double> 
	    S1ab(S1[sa_idx][fmt_.index_r(ej[0],sb,site)],
		 S1[sa_idx][fmt_.index_i(ej[0],sb,site)]);
	  */
	  complex<double> 
	    S1ab(S1[sa_idx][base_idx+ 2*Nc_*sb],
		 S1[sa_idx][base_idx+ 2*Nc_*sb+1]);


	  complex<double> S3S2ba(0.0,0.0);
	  complex<double> S2S3ba(0.0,0.0);
	  GammaResult Gb = G(sb);
	  
#pragma unroll(ND_)
	  for(int sc=0; sc<Nd_; ++sc){
	    GammaResult Gc = G(sc);
	    //common factor that we do not want to recalculate for every site...
	    //fac is an Nd*Nd matrix
	    //instead of passing G pass the precalculated matrix
	    //complex<double> fac(Gb.facr*Gc.facr -Gb.faci*Gc.faci,
	    //			Gb.facr*Gc.faci +Gb.faci*Gc.facr);

	   
	    /*!@brief CgS3Cg should be transversed for spinor idx. */
	    complex<double> 
	      CgS3Cg(S3[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_r(eps_[j*3+1],Gb.spn,site)],
		     S3[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_i(eps_[j*3+1],Gb.spn,site)]);

	    complex<double> 
	      CgS2Cg(S2[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_r(eps_[j*3+1],Gb.spn,site)],
		     S2[Nc_*Gc.spn+eps_[i*3+1]][fmt_.index_i(eps_[j*3+1],Gb.spn,site)]);

	    CgS2Cg *= fac_matrix[sc+Nd_*sb];	    
	    CgS3Cg *= fac_matrix[sc+Nd_*sb];

	    
	    S3S2ba += CgS3Cg
	      *complex<double>(S2[Nc_*sc+eps_[i*3+2]][fmt_.index_r(eps_[j*3+2],sa,site)],
			       S2[Nc_*sc+eps_[i*3+2]][fmt_.index_i(eps_[j*3+2],sa,site)]);
	    S2S3ba -= CgS2Cg
	      *complex<double>(S3[Nc_*sc+eps_[i*3+2]][fmt_.index_r(eps_[j*3+2],sa,site)],
			       S3[Nc_*sc+eps_[i*3+2]][fmt_.index_i(eps_[j*3+2],sa,site)]);

	    

	  }
	  Tr1 += S1ab*S3S2ba;
	  Tr2 += S1ab*S2S3ba;
	}
      }
      correl += eps*r1*Tr1;
      correl += eps*r2*Tr2;

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

  int slice3d = SiteIndex::instance()->slsize(0,TDIR);

#pragma omp parallel for
  for (int t = 0; t< Nt_; t++) {
    for(int site=0; site<slice3d; ++site){
      int site3d = SiteIndex::instance()->slice_t(t,site);
      
      correl_local[t] += spTrTr(S1,S1,S2,G,site3d,ud);
      correl_local[t] +=   spTr(S1,S1,S2,G,site3d,ud);
    }
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

   int slice3d = SiteIndex::instance()->slsize(0,TDIR);

  vector<complex<double> > fac_matrix;

  for (int a=0; a< Nd_; a++){
    GammaResult Ga = G(a);
    for (int b=0; b< Nd_; b++){
      GammaResult Gb = G(b);
      fac_matrix.push_back(complex<double> (Ga.facr*Gb.facr -Ga.faci*Gb.faci,
					    Ga.facr*Gb.faci +Ga.faci*Gb.facr));
    }
  }






#pragma omp parallel for
   for (int t = 0; t< Nt_; t++) {
     for(int site=0; site<slice3d; ++site){
       int site3d = SiteIndex::instance()->slice_t(t,site);
       //correl_local[t] += r2*spTrTr(Sl_,Sl_,Sh_,G,site3d,ud);
       //correl_local[t] += r4*spTrTr(Sh_,Sl_,Sl_,G,site3d,ud);
       correl_local[t] += spTrTrOpt(Sl_,Sl_,Sh_,G,r2,r4,site3d,ud);
       correl_local[t] += spTrOpt(  Sl_,Sh_,Sl_,G,r4,r2,fac_matrix,site3d,ud);
       correl_local[t] += r4*spTr(  Sh_,Sl_,Sl_,G,site3d,ud);
       //correl_local[t] -= r2*spTr(  Sl_,Sl_,Sh_,G,site3d,ud);
     }
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

  int slice3d = SiteIndex::instance()->slsize(0,TDIR);

#pragma omp parallel for
  for (int t = 0; t< Nt_; t++) {
    for(int site=0; site<slice3d; ++site){
      int site3d = SiteIndex::instance()->slice_t(t,site);
      correl_local[t] +=  spTrTr(S,S,S,*Cg,site3d,ud);
      correl_local[t] += r2*spTr(S,S,S,*Cg,site3d,ud);
    }
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
  complex<double> r1(1.0);
  complex<double> r2(2.0);


  vector<complex<double> > fac_matrix;

  for (int a=0; a< Nd_; a++){
    GammaResult Ga = (*Cg)(a);
    for (int b=0; b< Nd_; b++){
      GammaResult Gb = (*Cg)(b);
      fac_matrix.push_back(complex<double> (Ga.facr*Gb.facr -Ga.faci*Gb.faci,
					    Ga.facr*Gb.faci +Ga.faci*Gb.facr));
    }
  }




  int slice3d = SiteIndex::instance()->slsize(0,TDIR);

#pragma omp parallel for
  for (int t = 0; t< Nt_; t++) {
    for(int site=0; site<slice3d; ++site){
      int site3d = SiteIndex::instance()->slice_t(t,site);
      
      correl_local[t] +=    spTrTrOpt(S2,S1,S1,*Cg,r1,r2,site3d,ud);
      //correl_local[t] += r2*spTrTr(S1,S1,S2,*Cg,site3d,ud);
      correl_local[t] += spTrOpt(S1,S1,S2,*Cg,r2,r2, fac_matrix,site3d,ud);
      //correl_local[t] -= r2*spTrOpt(  S1,S2,S1,*Cg,fac_matrix,site3d,ud);
    }
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

