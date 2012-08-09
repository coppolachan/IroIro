//----------------------------------------------------------------------
/*! @file  dirac_staggered.cpp
 *  @brief Implementation of the staggered-Dirac operator
 */
//----------------------------------------------------------------------
#include "dirac_staggered.hpp"
#include "Main/Geometry/mapping.hpp"
#include "Tools/sunMat.hpp"
#include "Tools/sunVec.hpp"
#include <iostream>
#include <valarray>

void Dirac_staggered::set_ksp(){
  for(int site=0; site<Nvol_; ++site){
    int gs = SiteIndex::get_gsite((this->*get_site)(site)); 
    int x = SiteIndex::instance()->g_x(gs);
    int y = SiteIndex::instance()->g_y(gs);
    int z = SiteIndex::instance()->g_z(gs);
    
    ksp_[Nvol_*YDIR +site] *= double(1-2*(x%2));
    ksp_[Nvol_*ZDIR +site] *= double(1-2*((x+y)%2));
    ksp_[Nvol_*TDIR +site] *= double(1-2*((x+y+z)%2));
  }
}

void Dirac_staggered::set_ustg(){
  for(int mu=0; mu<Ndim_; ++mu)
    for(int site=0; site<Nvol_; ++site)
      ustg_.set(gf_.islice(site,mu),
		ksp_[Nvol_*mu+site]
		*SUNmat(u_[gf_.islice((this->*get_site)(site),mu)])); 
}

void Dirac_staggered::mult_offdiag(Field w,const Field& f)const{
  w = 0.0;
  for(int mu=0; mu<Ndim_; ++mu){
    Field ft = (this->*shift_fw)(f,mu);

    for(int site=0; site<Nvol_; ++site) 
      w.add(ff_.islice(site),
	    SUNmat(ustg_[gf_.islice((this->*get_site)(site),mu)])
	    *SUNvec(ft[ff_.islice(site)]));
    
    for(int site=0; site<Nvol_; ++site) 
      ft.set(ff_.islice(site),
	     SUNmat(ustg_[gf_.islice((this->*get_site)(site),mu)]).dag()
	     *SUNvec(f[ff_.islice(site)]));
  
    w -= (this->*shift_bk)(ft,mu);     
  }
  w *= 0.5/mq_;
}

void Dirac_staggered::mult_full(Field w,const Field& f)const{
  mult_offdiag(w,f);
  w += f;
}
void Dirac_staggered::mult_full_dag(Field w,const Field& f)const{
  mult_offdiag(w,f);
  -w;
  w += f;
}

const Field Dirac_staggered::mult(const Field& f)const{
  Field w(f.size());
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_staggered::mult_dag(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core_dag)(w,f);
  return w;
}

const Field Dirac_staggered::
md_force(const Field& eta,const Field& zeta)const{
  
  Field fp(gf_.size());
  md_force_p(fp,eta,zeta);
  md_force_m(fp,eta,zeta);
  fp *= -mp_;
  return fp;
}

void Dirac_staggered::
md_force_p(Field& fce,const Field& eta,const Field& zeta)const{}

void Dirac_staggered::
md_force_m(Field& fce,const Field& eta,const Field& zeta)const{}
