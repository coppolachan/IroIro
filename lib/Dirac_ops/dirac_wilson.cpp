/*! @file dirac_wilson.cpp
 *  @brief Declaration of Dirac_Wilson class
 * Time-stamp: <2013-08-20 16:32:48 cossu>
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
  /*
  double* temp = (double*) malloc(ff_.size()*sizeof(double));
  gamma5core(w,f);
  BGWilson_Mult(temp, pU,w,-kpp_,BGWILSON_DIRAC);
  gamma5core(w,temp);
  free(temp);
  */
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

  long double timing;
  FINE_TIMING_START(timing);

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
      
      /*
      for(int site=is; site<is+ns; ++site){
	int gsite = (this->*gp)(site);
	int shift = 2*NC_*NC_*(gsite+Nvol_*mu);
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
	    fce_ptr[shift+ 2*(a*NC_+b)  ] += fre;
	    fce_ptr[shift+ 2*(a*NC_+b)+1] += fim;
	    //	    f.set(a,b,fre,fim);
	  }
	}
	//fce.add(gf_.cslice(0,gsite,mu),f.getva());
      } 
      */

    }
  }

  FINE_TIMING_END(timing);
  _Message(TIMING_VERB_LEVEL, "[Timing] - Dirac_Wilson::md_force_p"
           << " - total timing = "
           << timing << std::endl);

}

void Dirac_Wilson::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;

  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

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
}

const Field Dirac_Wilson::md_force(const Field& eta,const Field& zeta)const{
  Field fp(gf_.size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);
  fp *= -kpp_;
  return fp;
}

