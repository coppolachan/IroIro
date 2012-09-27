//--------------------------------------------------------------
/*!@file dirac_wilson_EvenOdd.cpp
  @brief Definition of Even Odd wilson operator
*/
//-------------------------------------------------------------
#include "dirac_staggered_EvenOdd.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include <iostream>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace Mapping;

void Dirac_staggered_EvenOdd::set_ksphase(){
  for(int hs=0; hs<Nvh_; ++hs){
    // initialization for the even sector
    {
      int gs = SiteIndex::instance()->get_gsite(esec(hs)); 
      int x = SiteIndex::instance()->g_x(gs);
      int y = SiteIndex::instance()->g_y(gs);
      int z = SiteIndex::instance()->g_z(gs);
    
      kse_[Nvh_*YDIR +hs] *= double(1-2*(x%2));
      kse_[Nvh_*ZDIR +hs] *= double(1-2*((x+y)%2));
      kse_[Nvh_*TDIR +hs] *= double(1-2*((x+y+z)%2));
    }
    // initialization for the odd sector
    {
      int gs = SiteIndex::instance()->get_gsite(osec(hs)); 
      int x = SiteIndex::instance()->g_x(gs);
      int y = SiteIndex::instance()->g_y(gs);
      int z = SiteIndex::instance()->g_z(gs);
    
      kso_[Nvh_*YDIR +hs] *= double(1-2*(x%2));
      kso_[Nvh_*ZDIR +hs] *= double(1-2*((x+y)%2));
      kso_[Nvh_*TDIR +hs] *= double(1-2*((x+y+z)%2));
    }
  }
}

void Dirac_staggered_EvenOdd::set_ustag(){
  for(int mu=0; mu<Ndim_; ++mu){
    for(int hs=0; hs<Nvh_; ++hs){
      ue_.data.set(gf_.islice(hs,mu),
		   kse_[Nvh_*mu+hs]*(*u_)[gff_.islice(esec(hs),mu)]); 
      uo_.data.set(gf_.islice(hs,mu),
		   kso_[Nvh_*mu+hs]*(*u_)[gff_.islice(osec(hs),mu)]); 
    }
  }
  bdry_->apply_bc(ue_,uo_);
}

#if 1
#include "dirac_staggered_EvenOdd_improved.code"
#endif 

#if 0
void Dirac_staggered_EvenOdd::multPeo(Field& we,const Field& fo,int mu)const{
  FermionField1sp ft = shiftField_eo(FermionField1sp(fo),mu,Forward());
  for(int hs=0; hs<Nvh_; ++hs) 
    we.add(ff_.islice(hs),(mat(ue_,hs,mu)*vec(ft,hs)).getva());
}

void Dirac_staggered_EvenOdd::multPoe(Field& wo,const Field& fe,int mu)const{
  FermionField1sp ft = shiftField_oe(FermionField1sp(fe),mu,Forward());
  for(int hs=0; hs<Nvh_; ++hs)
    wo.add(ff_.islice(hs),(mat(uo_,hs,mu)*vec(ft,hs)).getva());
}

const Field Dirac_staggered_EvenOdd::mult_eo(const Field& fo)const{
  Field we(fsize_);
  for(int mu=0; mu<Ndim_; ++mu){
    multPeo(we,fo,mu);           //forward differenciation
    
    FermionField1sp ft(Nvh_);    
    for(int hs=0; hs<Nvh_; ++hs) //backward differenciation
      ft.data.set(ff_.islice(hs),
		  (mat_dag(uo_,hs,mu)*SUNvec(fo[ff_.islice(hs)])).getva());
    we -= shiftField_eo(ft,mu,Backward()).data;
  }
  we *= 0.5/mq_;
  return we;
}

const Field Dirac_staggered_EvenOdd::mult_oe(const Field& fe)const{
  Field wo(fsize_);
  for(int mu=0; mu<Ndim_; ++mu){
    multPoe(wo,fe,mu);           //forward differenciation
    
    FermionField1sp ft(Nvh_);    
    for(int hs=0; hs<Nvh_; ++hs) //backward differenciation
      ft.data.set(ff_.islice(hs),
		  (mat_dag(ue_,hs,mu)*SUNvec(fe[ff_.islice(hs)])).getva());

    wo -= shiftField_oe(ft,mu,Backward()).data;
  }
  wo *= 0.5/mq_;
  return wo;
}
#endif 0

const Field Dirac_staggered_EvenOdd::mult_eo_dag(const Field& fe) const{
  Field wo = mult_oe(fe);
  return -wo;
}

const Field Dirac_staggered_EvenOdd::mult_oe_dag(const Field& fo) const{
  Field we = mult_eo(fo);
  return -we;
}

const Field Dirac_staggered_EvenOdd::mult(const Field& fe) const{
  Field we(fe);
  we -= mult_eo(mult_oe(fe));
  return we;
}

const Field Dirac_staggered_EvenOdd::mult_dag(const Field& fe) const{
  return mult(fe);
}

const Field Dirac_staggered_EvenOdd::mult_full(const Field& f) const{
  assert(f.size()==2*fsize_);

  Field w(f);
  w.add(sle_,mult_eo(Field(f[slo_])).getva());
  w.add(slo_,mult_oe(Field(f[sle_])).getva());
  return w;
}

const Field Dirac_staggered_EvenOdd::mult_full_dag(const Field& f) const{
  assert(f.size()==2*fsize_);
  
  Field w(f);
  w.add(sle_,-(mult_eo(Field(f[slo_]))).getva());
  w.add(slo_,-(mult_oe(Field(f[sle_]))).getva());
  return w;
}

const Field Dirac_staggered_EvenOdd::
md_force(const Field& eta,const Field& zeta) const{
  Field fce(gsize_);

  for(int mu=0; mu<Ndim_; ++mu){
    Field etah(fsize_),zetah(fsize_);
    multPoe( etah, eta,mu);
    multPeo(zetah,zeta,mu);

    for(int hs=0; hs<Nvh_; ++hs){
      std::slice xsl = ff_.islice(hs);
      fce.set(gff_.islice(esec(hs),mu),
	      outer_prod_t(SUNvec(eta[xsl]),SUNvec(zetah[xsl])).getva());
      fce.set(gff_.islice(osec(hs),mu),
	      outer_prod_t(SUNvec(etah[xsl]),SUNvec(zeta[xsl])).getva());
    } 
  }
  fce *= -0.5/mq_;
  return fce;
}
