/*! @file dirac_wilson.cpp
 *  @brief Declaration of Dirac_Wilson class
 * Time-stamp: <2013-10-17 13:33:15 cossu>
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
    gamma5core(w.getaddr(ff_.index(0,site)),
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

#endif

/*!
 *  @brief MD-force contribution: \f$\zeta^\dagger\frac{dH_W}{d\tau}\eta\f$
 */
void Dirac_Wilson::md_force_p(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;
 

  double* fce_ptr  = fce.getaddr(0);
  double* eta_ptr  = const_cast<Field&>(eta).getaddr(0);
  double* zeta_ptr = const_cast<Field&>(zeta).getaddr(0);
  double* pU = const_cast<Field *>(u_)->getaddr(0);
  long double timing;
  FINE_TIMING_START(timing);
#ifdef IBM_BGQ_WILSON
  Field xie(ff_.size());
  double* xie_ptr  = xie.getaddr(0);
  vector_int global_sites;
  if (EO_BGWilson == 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else
    global_sites = SiteIndex_EvenOdd::instance()->osec();

#pragma omp parallel 
  { 
    //SUNmat f;
    int tid, nid;
    int is, ie, ns;
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = tid*Nvol_ / nid;
    ns = Nvol_ / nid;
    
    for(int mu=0; mu<NDIM_; ++mu){
      BGWilson_MultEO_Dir(xie_ptr, pU, eta_ptr,  1.0, EO_BGWilson, BGWILSON_DIRAC, mu, BGWILSON_FORWARD);
      
      for(int site=is; site<is+ns; ++site){
	//f = 0.0;
	unsigned int index = ff_.Nin()*site;
	unsigned int g_idx = gf_.index(0,global_sites[site],mu); 
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    unsigned int fce_idx = g_idx+2*(NC_*a+b);
	    for(int s=0; s<ND_; ++s){
	      unsigned int ra = index+ 2*(a+NC_*s);
	      unsigned int rb = index+ 2*(b+NC_*s);
	   	      
	      fce_ptr[fce_idx  ] += zeta[rb]*xie[ra  ] +zeta[rb+1]*xie[ra+1];
	      fce_ptr[fce_idx+1] += zeta[rb]*xie[ra+1] -zeta[rb+1]*xie[ra  ];
	    }//spin
	  }//b
	}//a 
      } //site
    }//mu
  }//omp
  
  FINE_TIMING_END(timing);
  _Message(DEBUG_VERB_LEVEL, "[Timing] - Dirac_Wilson::md_force_p"
           << " - total timing = "
           << timing << std::endl);
#else
  for(int mu=0; mu<NDIM_; ++mu){
    Field xie(ff_.size());
    (this->*mult_p[mu])(xie, eta);
       
#pragma omp parallel 
    {
      int tid, nid;
      int is, ie, ns;
      nid = omp_get_num_threads();
      tid = omp_get_thread_num();
      is = tid*Nvol_ / nid;
      ns = Nvol_ / nid;
      
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

  FINE_TIMING_END(timing);
  _Message(DEBUG_VERB_LEVEL, "[Timing] - Dirac_Wilson::md_force_p"
           << " - total timing = "
           << timing << std::endl);
#endif

}




void Dirac_Wilson::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;

  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  double* fce_ptr  = fce.getaddr(0);
  double* et5_ptr  = et5.getaddr(0);
  double* zt5_ptr  = zt5.getaddr(0);
  double* pU = const_cast<Field *>(u_)->getaddr(0);

#ifdef IBM_BGQ_WILSON
  Field xz5(ff_.size());
  double* xz5_ptr  = xz5.getaddr(0);
  vector_int global_sites;
  if (EO_BGWilson == 1)
    global_sites = SiteIndex_EvenOdd::instance()->esec();
  else
    global_sites = SiteIndex_EvenOdd::instance()->osec();

#pragma omp parallel 
  { 
    //SUNmat f;
    int tid, nid;
    int is, ie, ns;
    nid = omp_get_num_threads();
    tid = omp_get_thread_num();
    is = tid*Nvol_ / nid;
    ns = Nvol_ / nid;
    
    for(int mu=0; mu<NDIM_; ++mu){
      BGWilson_MultEO_Dir(xz5_ptr, pU, zt5_ptr,  1.0, EO_BGWilson, BGWILSON_DIRAC, mu, BGWILSON_FORWARD);
      
      for(int site=is; site<is+ns; ++site){
	//f = 0.0;
	unsigned int index = ff_.Nin()*site;
	unsigned int g_idx = gf_.index(0,global_sites[site],mu); 
	for(int a=0; a<NC_; ++a){
	  for(int b=0; b<NC_; ++b){
	    double fre = 0.0;
	    double fim = 0.0;

	    for(int s=0; s<ND_; ++s){
	      unsigned int ra = index+ 2*(a+NC_*s);
	      unsigned int rb = index+ 2*(b+NC_*s);
	   	
	      fre -= xz5[rb]*et5[ra  ] +xz5[rb+1]*et5[ra+1];
	      fim -= xz5[rb]*et5[ra+1] -xz5[rb+1]*et5[ra  ];
	    }//spin
	    fce_ptr[g_idx+2*(NC_*a+b)  ] += fre;
	    fce_ptr[g_idx+2*(NC_*a+b)+1] += fim;
	  }//b
	}//a 
      } //site
      
    }//mu
  }//omp

#else
  for(int mu=0; mu<NDIM_; ++mu){
    Field xz5(ff_.size());
    (this->*mult_p[mu])(xz5, zt5);
    

#pragma omp parallel 
    {
      int tid, nid;
      int is, ie, ns;
      nid = omp_get_num_threads();
      tid = omp_get_thread_num();
      is = tid*Nvol_ / nid;
      ie = (tid + 1)*Nvol_ / nid;
      ns = ie - is;
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

