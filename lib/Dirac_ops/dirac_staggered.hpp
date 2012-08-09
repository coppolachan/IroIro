//-------------------------------------------------------------
/*!@file   dirac_staggered.hpp
 * @brief  Definition of the staggered operator
 */
//-------------------------------------------------------------
#ifndef DIRAC_STAGGERED_INCLUDED
#define DIRAC_STAGGERED_INCLUDED

#include "dirac.hpp"
#include "include/format_S.h"
#include "include/format_G.h"
#include "Main/Geometry/mapping.hpp"

typedef Format::Format_S sfmt_t;
typedef Format::Format_G gfmt_t;

class Dirac_staggered: public DiracStaggeredLike{

private:
  const Field* const u_;
  Field ustg_;
  int Nvol_,Ndim_;
  double mq_;
  std::vector<int> bdry_;
  const size_t fsize_;
  const size_t gsize_;
  std::valarray<double> ksp_; /*!@brief Kawamoto-Schmidt phase*/

  const sfmt_t ff_;
  const gfmt_t gf_;

  void set_ksp();
  void set_ustg();

  //getting site idx
  int(Dirac_staggered::*get_site)(int)const;
  int get_fullsite(int site)const{return site;}
  int get_esite(int hs)const{return SiteIndex_EvenOdd::instance()->esec(hs);}
  int get_osite(int hs)const{return SiteIndex_EvenOdd::instance()->osec(hs);}

  //field shift (forward)
  const Field(Dirac_staggered::*shift_fw)(const Field&,int)const;
  const Field shift_full_fw(const Field& f,int mu)const{
    return Mapping::shiftField(f,mu,Forward());
  }
  const Field shift_oe_fw(const Field& f,int mu)const{
    return Mapping::shiftField_oe(f,mu,Forward());
  }
  const Field shift_eo_fw(const Field& f,int mu)const{
    return Mapping::shiftField_eo(f,mu,Forward());
  }
  
  //field shift (backward)
  const Field(Dirac_staggered::*shift_bk)(const Field&,int)const;
  const Field shift_full_bk(const Field& f,int mu)const{
    return Mapping::shiftField(f,mu,Backrward());
  }
  const Field shift_oe_bk(const Field& f,int mu)const{
    return Mapping::shiftField_oe(f,mu,Backrward());
  }
  const Field shift_eo_bk(const Field& f,int mu)const{
    return Mapping::shiftField_eo(f,mu,Backrward());
  }

  //main part of mult 
  void(Dirac_staggered::*mult_core)(Field&,const Field&)const;
  void mult_full(   Field&,const Field&)const;  /*! @brief  (1-kpp*D)*f */
  void mult_offdiag(Field&,const Field&)const;  /*! @brief  -kpp*D*f */

  void(Dirac_staggered::*mult_core_dag)(Field&,const Field&)const;

public:
  // full-site-index (manual instanciation)
  Dirac_staggered(double mass,const Field* u)   
    :mq_(mass),u_(u),ustg_(u->size()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(Ndim_),
     ff_(Nvol_),gf_(Nvol_),
     ksp_(1.0,Nvol_*Ndim_),
     fsize_(ff_.size()),gsize_(gf_.size()),
     get_site(&Dirac_staggered::get_fullsite),
     shift_fw(&Dirac_staggered::shift_full_fw),
     shift_bk(&Dirac_staggered::shift_full_bk),
     mult_core(&Dirac_staggered::mult_full),
     mult_core_dag(&Dirac_staggered::mult_full){
    //
    set_ksp();
    set_ustg();
  }

  // full-site-index (xml instanciation)
  Dirac_staggered(const XML::node& node,const Field* u)
    :u_(u),ustg_(u->size()),
     Nvol_(CommonPrms::instance()->Nvol()),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(Ndim_),
     ff_(Nvol_),gf_(Nvol_),
     ksp_(1.0,Nvol_*Ndim_),
     fsize_(ff_.size()),gsize_(gf_.size()),
     get_site(&Dirac_staggered::get_fullsite),
     shift_fw(&Dirac_staggered::shift_full_fw),
     shift_bk(&Dirac_staggered::shift_full_bk),
     mult_core(&Dirac_staggered::mult_full),
     mult_core_dag(&Dirac_staggered::mult_full){
    //
    XML::read(node,"mass",mq_);
    set_ksp();
    set_ustg();
  }

  // odd-to-even (Deo) 
  Dirac_staggered(double mass,const Field* u,Dop::EOtag)   
    :mq_(mass),u_(u),ustg_(u->size()),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(Ndim_),
     ff_(Nvol_),gf_(Nvol_),
     ksp_(1.0,Nvol_*Ndim_),
     fsize_(ff_.size()),gsize_(gf_.size()),
     get_site(&Dirac_staggered::get_osite),
     shift_fw(&Dirac_staggered::shift_eo_fw),
     shift_bk(&Dirac_staggered::shift_eo_bk),
     mult_core(&Dirac_staggered::mult_offdiag),
     mult_core_dag(&Dirac_staggered::mult_offdiag){
    //
    set_ksp();
    set_ustg();
  }

  // even-to-odd (Doe)
  Dirac_staggered(double mass,const Field* u,Dop::OEtag)   
    :mq_(mass),u_(u),ustg_(u->size()),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(Ndim_),
     ff_(Nvol_),gf_(Nvol_),
     ksp_(1.0,Nvol_*Ndim_),
     fsize_(ff_.size()),gsize_(gf_.size()),
     get_site(&Dirac_staggered::get_esite),
     shift_fw(&Dirac_staggered::shift_oe_fw),
     shift_bk(&Dirac_staggered::shift_oe_bk),
     mult_core(&Dirac_staggered::mult_offdiag),
     mult_core_dag(&Dirac_staggered::mult_offdiag){
    //
    set_ksp();
    set_ustg();
  }

  ~Dirac_staggered(){}

  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // EvenOdd operator has no preconditioner now 
  const Field mult_prec     (const Field& f)const{return f;}
  const Field mult_dag_prec (const Field& f)const{return f;}
  const Field left_prec     (const Field& f)const{return f;}
  const Field right_prec    (const Field& f)const{return f;}
  const Field left_dag_prec (const Field& f)const{return f;}
  const Field right_dag_prec(const Field& f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field md_force(const Field&,const Field&) const;
  void md_force_p(Field&,const Field&,const Field&)const;
  void md_force_m(Field&,const Field&,const Field&)const;

  const std::vector<int> get_gsite() const{
    return SiteIndex::instance()->get_gsite();
  }
  
  double get_mq()const{ return mq_;}

  void update_internal_state(){ set_ustg();}
};

#endif
