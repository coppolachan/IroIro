/*! @file dirac_wilson.cpp
 *  @brief Declaration of Dirac_Wilson class
 * Time-stamp: <2014-08-09 14:28:41 noaki>
 */
#include "dirac_wilson.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include <omp.h>

#ifdef IBM_BGQ_WILSON
#include "bgqwilson.h"
#include "bgqthread.h"
#include "Tools/Architecture_Optimized/utils_BGQ.hpp"
#include "Architecture_Optimized/dirac_wilson_BGQ.code"

#elif defined IMPROVED_WILSON
#include "dirac_wilson_improved.code"
#else
#include "dirac_wilson_standard.code"
#endif 

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

using namespace SUNvecUtils;
using namespace std;

const Field Dirac_Wilson::mult(const Field& f) const{
  Field w(ff_.size());
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_Wilson::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

/*!
 *  @brief MD-force contribution: \f$\zeta^\dagger\frac{dH_W}{d\tau}\eta\f$
 */
void Dirac_Wilson::md_force_p(Field& fce,const Field& eta,const Field& zeta)const{
  for(int mu=0; mu<NDIM_; ++mu) mkfrc(fce,eta,zeta,mu);
}

void Dirac_Wilson::md_force_m(Field& fce,const Field& eta,const Field& zeta)const{
  Field zt5 = gamma5(zeta);
  Field et5 = gamma5(eta);
  for(int mu=0; mu<NDIM_; ++mu) mkfrc(fce,zt5,et5,mu);
}

const Field Dirac_Wilson::md_force(const Field& eta,const Field& zeta)const{
  Field fp(gf_.size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);

  fp *= -kpp_;
  return fp;
}

void Dirac_Wilson::mkfrc(Field& fce,const Field& eta,const Field& zeta,int mu)const{
  using namespace SUNmatUtils;

  Field xie(ff_.size());
  (this->*mult_p[mu])(xie,eta);
#pragma omp parallel 
  {
    int ns = Nvol_/omp_get_num_threads();
    int is = omp_get_thread_num()*ns;

    for(int site=is; site<is+ns; ++site){
      for(int a=0; a<NC_; ++a){
	for(int b=0; b<NC_; ++b){

	  for(int s=0; s<ND_; ++s){
	    size_t ra =ff_.index_r(a,s,site);
	    size_t ia =ff_.index_i(a,s,site);
	    
	    size_t rb =ff_.index_r(b,s,site);
	    size_t ib =ff_.index_i(b,s,site);

	    fce.add(gf_.index_r(a,b,gp_[site],mu), 
		    zeta[rb]*xie[ra]+zeta[ib]*xie[ia]);

	    fce.add(gf_.index_i(a,b,gp_[site],mu), 
		    zeta[rb]*xie[ia]-zeta[ib]*xie[ra]);
	  }
	}
      }
    }
  }
}

#ifndef IBM_BGQ_WILSON
//////////////// regular implementations ////////////////
void Dirac_Wilson::init_mult_pm(){
  mult_p[XDIR] = &Dirac_Wilson::mult_xp; 
  mult_p[YDIR] = &Dirac_Wilson::mult_yp; 
  mult_p[ZDIR] = &Dirac_Wilson::mult_zp; 
  mult_p[TDIR] = &Dirac_Wilson::mult_tp; 
	       
  mult_m[XDIR] = &Dirac_Wilson::mult_xm; 
  mult_m[YDIR] = &Dirac_Wilson::mult_ym; 
  mult_m[ZDIR] = &Dirac_Wilson::mult_zm; 
  mult_m[TDIR] = &Dirac_Wilson::mult_tm; 
}

void Dirac_Wilson::mult_offdiag(Field& w, const Field& f) const{
  for(int d=0; d <NDIM_; ++d){
    (this->*mult_p[d])(w,f);
    (this->*mult_m[d])(w,f);
  }
  w *= -kpp_;
}

void Dirac_Wilson::mult_full(Field& w, const Field& f) const{
  mult_offdiag(w,f);
  w += f; 
}

const Field Dirac_Wilson::gamma5(const Field& f)const{ 
  Field w(ff_.size());
  for(int site=0; site<Nvol_; ++site)
    dm_.gamma5core(w.getaddr(ff_.index(0,site)),
		   const_cast<Field&>(f).getaddr(ff_.index(0,site)));
  return w;
}

#endif

