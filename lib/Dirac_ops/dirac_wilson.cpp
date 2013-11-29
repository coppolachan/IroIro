/*! @file dirac_wilson.cpp
 *  @brief Declaration of Dirac_Wilson class
 * Time-stamp: <2013-11-27 17:20:38 noaki>
 */
#include "dirac_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"

using namespace SUNvecUtils;
using namespace std;

#include <omp.h>
#ifdef IMPROVED_WILSON
#include "dirac_wilson_improved.code"
#else
#include "dirac_wilson_standard.code"
#endif /*IMPROVED_WILSON*/

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#include "bgqthread.h"
#include "Tools/utils_BGQ.hpp"
static Field xie;
#endif

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

void (Dirac_Wilson::*Dirac_Wilson::mult_p[])
(Field&,const Field&) const = {&Dirac_Wilson::mult_xp,
			       &Dirac_Wilson::mult_yp,
			       &Dirac_Wilson::mult_zp,
			       &Dirac_Wilson::mult_tp,};

void (Dirac_Wilson::*Dirac_Wilson::mult_m[])
(Field&,const Field&) const = {&Dirac_Wilson::mult_xm,
			       &Dirac_Wilson::mult_ym,
			       &Dirac_Wilson::mult_zm,
			       &Dirac_Wilson::mult_tm,};

void Dirac_Wilson::mult_offdiag(Field& w, const Field& f) const{
#ifndef IBM_BGQ_WILSON
  for(int d=0; d <NDIM_; ++d){
    (this->*mult_p[d])(w,f);
    (this->*mult_m[d])(w,f);
  }
  w *= -kpp_;
#else  
  double* pF = const_cast<Field&>(f).getaddr(0);
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  double* pW = w.getaddr(0);
  BGWilson_MultEO(pW, pU, pF, -kpp_ , EO_BGWilson, BGWILSON_DIRAC);
#endif
}
void Dirac_Wilson::mult_full(Field& w, const Field& f) const{
#ifndef IBM_BGQ_WILSON
  mult_offdiag(w,f);
  w += f; 
#else  
  double* pF = const_cast<Field&>(f).getaddr(0);
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  double* pW = w.getaddr(0);
  BGWilson_Mult(pW, pU, pF, -kpp_ , BGWILSON_DIRAC);
#endif
}

const Field Dirac_Wilson::mult(const Field& f) const{
  Field w(ff_.size());
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_Wilson::gamma5(const Field& f)const{ 
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site){
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  }
  return w;
}

const Field Dirac_Wilson::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

#ifdef IBM_BGQ_WILSON
void Dirac_Wilson::mult_ptr(double* w, double* const f) const{
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  BGWilson_Mult(w,pU,f,-kpp_,BGWILSON_DIRAC);
}
void Dirac_Wilson::mult_dag_ptr(double* w, double* const f) const{
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  BGWilson_Mult_Dag(w,pU,f,-kpp_,BGWILSON_DIRAC);
}
void Dirac_Wilson::mult_ptr_EO(double* w, double* const f) const{
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  BGWilson_MultEO(w,pU,f,-kpp_,EO_BGWilson,BGWILSON_DIRAC);
}
void Dirac_Wilson::mult_dag_ptr_EO(double* w, double* const f) const{
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  BGWilson_MultEO_Dag(w,pU,f,-kpp_,EO_BGWilson,BGWILSON_DIRAC);
}
void Dirac_Wilson::BGQ_initialize_pointers(){
  xie.resize(ff_.size());//(double*)malloc(ff_.Nin()*Nvol_*sizeof(double));
  if (EO_BGWilson == 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else
    global_sites = SiteIndex_EvenOdd::instance()->osec();
}
#endif

/*!
 *  @brief MD-force contribution: \f$\zeta^\dagger\frac{dH_W}{d\tau}\eta\f$
 */
void Dirac_Wilson::md_force_p(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;
 
#ifdef IBM_BGQ_WILSON
  double* fce_ptr  = fce.getaddr(0);
  double* eta_ptr  = const_cast<Field&>(eta).getaddr(0);
  double* zeta_ptr = const_cast<Field&>(zeta).getaddr(0);
  double* pU       = const_cast<Field *>(u_)->getaddr(0);
  double* xie_ptr  = const_cast<Field&>(xie).getaddr(0);

  //#pragma omp parallel 
  //{
    int nid = omp_get_num_threads();
    int ns = Nvol_/nid;
    int is = omp_get_thread_num()*ns;

    for(int mu=0; mu<NDIM_; ++mu){
      BGWilson_MultEO_Dir(xie_ptr, pU, eta_ptr, 1.0, EO_BGWilson, BGWILSON_DIRAC, mu, BGWILSON_FORWARD);
      
      for(int site=is; site<is+ns; ++site){
	unsigned int index = ff_.Nin()*site;
	unsigned int g_idx = gf_.index(0,global_sites[site],mu); 
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    unsigned int fce_idx = g_idx+2*(NC_*a+b);
	    for(int s=0; s<ND_; ++s){
	      unsigned int ra = index+ 2*(a+NC_*s);
	      unsigned int rb = index+ 2*(b+NC_*s);
	      
	      fce_ptr[fce_idx  ] += zeta_ptr[rb]*xie_ptr[ra  ] +zeta_ptr[rb+1]*xie_ptr[ra+1];
	      fce_ptr[fce_idx+1] += zeta_ptr[rb]*xie_ptr[ra+1] -zeta_ptr[rb+1]*xie_ptr[ra  ];
	    }//spin
	  }//b
	}//a 
      } //site
    }//mu
     
    //BGQThread_Barrier(0,nid);
    //}//omp
  
#else
  for(int mu=0; mu<NDIM_; ++mu){
    Field xie(ff_.size());
    (this->*mult_p[mu])(xie, eta);
       
#pragma omp parallel 
    {
      int ns = Nvol_/omp_get_num_threads();
      int is = omp_get_thread_num()*ns;
      SUNmat f;
      
      for(int site=is; site<is+ns; ++site){
	f = 0.0;
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    double fre = 0.0;
	    double fim = 0.0;
	    for(int s=0; s<ND_; ++s){
	      
	      size_t ra =ff_.index_r(a,s,site);
	      size_t ia =ff_.index_i(a,s,site);
	      
	      size_t rb =ff_.index_r(b,s,site);
	      size_t ib =ff_.index_i(b,s,site);
	      
	      fre += zeta[rb]*xie[ra] +zeta[ib]*xie[ia];
	      fim += zeta[rb]*xie[ia] -zeta[ib]*xie[ra];
	    }
	    f.set(a,b,fre,fim);
	  }
	}
	int gsite = (this->*gp)(site);
	fce.add(gf_.cslice(0,gsite,mu),f.getva());
      } 
    }
  }
#endif
}

void Dirac_Wilson::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;

  Field zt5(ff_.size());
  Field et5(ff_.size());

#ifdef IBM_BGQ_WILSON
  double* fce_ptr  = fce.getaddr(0);
  double* et5_ptr  = et5.getaddr(0);
  double* zt5_ptr  = zt5.getaddr(0);
  Spinor* zeta_ptr = (Spinor*)const_cast<Field&>(zeta).getaddr(0);
  Spinor* eta_ptr  = (Spinor*)const_cast<Field&>(eta).getaddr(0);
  double* pU = const_cast<Field *>(u_)->getaddr(0);

  Field xz5(ff_.size());
  double* xz5_ptr  = xz5.getaddr(0);
  vector_int global_sites;
  if (EO_BGWilson == 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else
    global_sites = SiteIndex_EvenOdd::instance()->osec();

#pragma omp parallel 
  { 
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;
    
    BGWilsonLA_MultGamma5((Spinor*)(zt5_ptr)+is, zeta_ptr+is, ns);
    BGWilsonLA_MultGamma5((Spinor*)(et5_ptr)+is, eta_ptr+is, ns);

    for(int mu=0; mu<NDIM_; ++mu){
      BGWilson_MultEO_Dir(xz5_ptr, pU, zt5_ptr,  1.0, EO_BGWilson, BGWILSON_DIRAC, mu, BGWILSON_FORWARD);
      
      for(int site=is; site<is+ns; ++site){
	//f = 0.0;
	unsigned int index = ff_.Nin()*site;
	unsigned int g_idx = gf_.index(0,global_sites[site],mu); 
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    //double fre = 0.0;
	    //double fim = 0.0;
	    unsigned int fce_idx = g_idx+2*(NC_*a+b);
	    for(int s=0; s<ND_; ++s){
	      unsigned int ra = index+ 2*(a+NC_*s);
	      unsigned int rb = index+ 2*(b+NC_*s);
	   	
	      fce_ptr[fce_idx  ] -= xz5[rb]*et5[ra  ] +xz5[rb+1]*et5[ra+1];
	      fce_ptr[fce_idx+1] -= xz5[rb]*et5[ra+1] -xz5[rb+1]*et5[ra  ];
	    }//spin
	    //fce_ptr[g_idx+2*(NC_*a+b)  ] += fre;
	    //fce_ptr[g_idx+2*(NC_*a+b)+1] += fim;
	  }//b
	}//a 
      } //site
      
    }//mu
  }//omp

#else
  zt5 = gamma5(zeta);
  et5 = gamma5(eta);

  for(int mu=0; mu<NDIM_; ++mu){
    Field xz5(ff_.size());
    (this->*mult_p[mu])(xz5, zt5);

#pragma omp parallel 
    {
      int ns = Nvol_/omp_get_num_threads();
      int is = omp_get_thread_num()*ns;
      SUNmat f;
      
      for(int site=is; site<is+ns; ++site){
	f=0.0;
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    double fre = 0.0;
	    double fim = 0.0;

	    for(int s=0; s<ND_; ++s){
	    
	      size_t ra =ff_.index_r(a,s,site);
	      size_t ia =ff_.index_i(a,s,site);
	    
	      size_t rb =ff_.index_r(b,s,site);
	      size_t ib =ff_.index_i(b,s,site);
	    
	      fre -= xz5[rb]*et5[ra] +xz5[ib]*et5[ia];
	      fim -= xz5[rb]*et5[ia] -xz5[ib]*et5[ra];
	    }
	    f.set(a,b,fre,fim);
	  }
	}
	int gsite = (this->*gp)(site);
	fce.add(gf_.cslice(0,gsite,mu),f.getva());
      }
    }
  }
#endif

}

const Field Dirac_Wilson::md_force(const Field& eta,const Field& zeta)const{
  Field fp(gf_.size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);
  fp *= -kpp_;
  return fp;
}

