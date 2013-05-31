//--------------------------------------------------------------
/*!@file dirac_wilson_EvenOdd_adjoint.cpp
  @brief Definition of the Even Odd staggered operator
*/
//-------------------------------------------------------------

#include "dirac_staggered_EvenOdd.hpp"
#include "Fields/field_expressions.hpp"
#include "Tools/fieldUtils.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/sunVec.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include <iostream>

using namespace std;
using namespace FieldUtils;
using namespace SUNmatUtils;
using namespace SUNvecUtils;
using namespace Mapping;

void (Dirac_staggered_EvenOdd::*Dirac_staggered_EvenOdd::mult_type[])
(Field&,const Field&) const = {&Dirac_staggered_EvenOdd::mult_DdagDee,
			       &Dirac_staggered_EvenOdd::mult_DdagDoo,
			       &Dirac_staggered_EvenOdd::mult_Dfull,};

void (Dirac_staggered_EvenOdd::*Dirac_staggered_EvenOdd::mult_dag_type[])
(Field&,const Field&) const = {&Dirac_staggered_EvenOdd::mult_DdagDee,
			       &Dirac_staggered_EvenOdd::mult_DdagDoo,
			       &Dirac_staggered_EvenOdd::mult_Dfull_dag,};

void Dirac_staggered_EvenOdd::set_ustag(){
  /*
  Staples stpl;
  double plaq = stpl.plaquette(GaugeField(*u_));
  CCIO::cout<<"Dirac_staggered_EvenOdd::set_ustag   plaq="<<plaq<<"\n";
  */
  for(int mu=0; mu<Ndim_; ++mu){
    for(int hs=0; hs<Nvh_; ++hs){
      ue_.data.set(gf_.islice(hs,mu),
		   kse_[Nvh_*mu+hs]
		   *(*u_)[gff_.islice(SiteIndex_EvenOdd::instance()->esec(hs),
				      mu)]); 
      uo_.data.set(gf_.islice(hs,mu),
		   kso_[Nvh_*mu+hs]
		   *(*u_)[gff_.islice(SiteIndex_EvenOdd::instance()->osec(hs),
				      mu)]); 
    }
  }
  ue_*= 0.5/mq_;
  uo_*= 0.5/mq_;
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

    double ft_norm = ft.data.norm();
    Field sft = shiftField_eo(ft,mu,Backward()).data;

    std::cout.precision(16);
    for(int os=0;os<Nvh_;++os){
      int idx = ff_.index(0,os);
      int site = SiteIndex::instance()->get_gsite(SiteIndex_EvenOdd::instance()->osec(os));
      int ix_o = SiteIndex::instance()->g_x(site);
      int iy_o = SiteIndex::instance()->g_y(site);
      int iz_o = SiteIndex::instance()->g_z(site);
      int it_o = SiteIndex::instance()->g_t(site);
      
      site = SiteIndex::instance()->get_gsite(SiteIndex_EvenOdd::instance()->esec(os));
      int ix_e = SiteIndex::instance()->g_x(site);
      int iy_e = SiteIndex::instance()->g_y(site);
      int iz_e = SiteIndex::instance()->g_z(site);
      int it_e = SiteIndex::instance()->g_t(site);
      /*
	std::cout<<ix<<","<<iy<<","<<iz<<","<<it
	  	 <<" sft["<<idx<<"]="<<sft[idx]
		 <<" we["<<idx<<"]="<<we[idx]<<"\n";
      */
      std::cout<<ix_o<<","<<iy_o<<","<<iz_o<<","<<it_o
	       <<" ft["<<idx<<"]="<<ft[idx]<<"  "
	       <<ix_e<<","<<iy_e<<","<<iz_e<<","<<it_e
	       <<" sft["<<idx<<"]="<<sft[idx]<<"\n";
    }

    double sft_norm = sft.norm();

    double w_norm0 = we.norm();
    we -= sft;

    //we -= shiftField_eo(ft,mu,Backward()).data;

    double w_norm = we.norm();

    CCIO::cout<<" ft_norm=" <<ft_norm
	      <<" sft_norm="<<sft_norm
	      <<" w_norm0="<< w_norm0
	      <<" w_norm="<< w_norm<<"\n";
  }
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
  return wo;
}
#endif //0

const Field Dirac_staggered_EvenOdd::mult_eo_dag(const Field& fe) const{
  return -Field(mult_oe(fe));
}

const Field Dirac_staggered_EvenOdd::mult_oe_dag(const Field& fo) const{
  return -Field(mult_eo(fo));
}

const Field Dirac_staggered_EvenOdd::mult(const Field& f)const{
  Field w(fsize_);
  (this->*mult_core)(w,f);
  return w;
}

const Field Dirac_staggered_EvenOdd::mult_dag(const Field& f)const{
  Field w(fsize_);
  (this->*mult_dag_core)(w,f);
  return w;
}

void Dirac_staggered_EvenOdd::mult_DdagDee(Field& we,const Field& fe)const{
  using namespace FieldExpression;
  we = fe -mult_eo(mult_oe(fe));
}

void Dirac_staggered_EvenOdd::mult_DdagDoo(Field& wo,const Field& fo)const{
  using namespace FieldExpression;
  wo = fo -mult_oe(mult_eo(fo));
}

void Dirac_staggered_EvenOdd::mult_Dfull(Field& w,const Field& f)const{
  assert(f.size()==2*fsize_);
  assert(w.size()==2*fsize_);
  w = f;
  w.add(sle_,mult_eo(Field(f[slo_])).getva());
  w.add(slo_,mult_oe(Field(f[sle_])).getva());
}

void Dirac_staggered_EvenOdd::mult_Dfull_dag(Field& w,const Field& f) const{
  assert(f.size()==2*fsize_);
  assert(w.size()==2*fsize_);
  w = f;
  w.add(sle_,-(mult_eo(Field(f[slo_]))).getva());
  w.add(slo_,-(mult_oe(Field(f[sle_]))).getva());
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
      fce.set(gff_.islice(SiteIndex_EvenOdd::instance()->esec(hs),mu),
	      outer_prod_t(SUNvec(eta[xsl]),SUNvec(zetah[xsl])).getva());
      fce.set(gff_.islice(SiteIndex_EvenOdd::instance()->osec(hs),mu),
	      outer_prod_t(SUNvec(etah[xsl]),SUNvec(zeta[xsl])).getva());
    } 
  }
  fce *= -1.0;
  return fce;
}


