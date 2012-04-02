//----------------------------------------------------------------------
// dirac_wilson.cpp
//----------------------------------------------------------------------
#include "dirac_wilson.hpp"
#include "Tools/sunMatUtils.hpp"

using namespace SUNvec_utils;
using namespace std;

#ifdef IMPROVED_WILSON
#include "dirac_wilson_improved.code"
#else
#include "dirac_wilson_standard.code"
#endif /*IMPROVED_WILSON*/

/////////////////////////////////////////////////////////////////////////

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

const Field Dirac_Wilson::mult(const Field& f) const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_Wilson::mult_dag(const Field& f)const{ 
  return gamma5(mult(gamma5(f)));
}

/*!
 *  @brief MD-force contribution: \f$\zeta^\dagger\frac{dH_W}{d\tau}\eta\f$
 */
void Dirac_Wilson::md_force_p(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;
  SUNmat f;

  for(int mu=0; mu<NDIM_; ++mu){
    Field xie(fsize_);

    (this->*mult_p[mu])(xie, eta);

    for(int site=0; site<Nvol_; ++site){
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

void Dirac_Wilson::md_force_m(Field& fce,
			      const Field& eta,const Field& zeta)const{
  using namespace SUNmatUtils;
  SUNmat f;
  Field et5 = gamma5(eta);
  Field zt5 = gamma5(zeta);

  for(int mu=0; mu<NDIM_; ++mu){
    Field xz5(fsize_);
    (this->*mult_p[mu])(xz5, zt5);

    for(int site=0; site<Nvol_; ++site){
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

const Field Dirac_Wilson::
md_force(const Field& eta,const Field& zeta)const{
  
  Field fp(gf_.size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);
  fp *= -kpp_;
  return fp;
}

const vector<int> Dirac_Wilson::get_gsite() const {
  return SiteIndex::instance()->get_gsite();
}

namespace Dw{
  double read_mass(const XML::node& node){
    double mass;
    XML::read(node, "mass", mass);
    return mass;
  }
}
