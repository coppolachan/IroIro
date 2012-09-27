//-------------------------------------------------------------
/*!@file   dirac_staggered_EvenOdd.hpp
 * @brief  Definition of the staggered operator
 */
//-------------------------------------------------------------
#ifndef DIRAC_STAGGERED_EVENODD_INCLUDED
#define DIRAC_STAGGERED_EVENODD_INCLUDED

#include "dirac.hpp"
#include "include/format_S.h"
#include "include/format_G.h"
#include "Main/Geometry/shiftField.hpp"
#include "boundaryCond.hpp"
#include <memory>

typedef Format::Format_S sfmt_t;
typedef Format::Format_G gfmt_t;

class Dirac_staggered_EvenOdd :public DiracStaggeredEvenOddLike{
private:
  const Field* const u_;
  int Nvh_,Ndim_;
  double mq_;
  GaugeField ue_,uo_;

  const sfmt_t ff_;  // for all pseudo-fermions 
  const gfmt_t gf_;  // for ustg_
  const gfmt_t gff_; // for u_ (always full site)

  const size_t fsize_;
  const size_t gsize_;

  const BoundaryCond* bdry_; 
  
  std::valarray<double> kse_;  /*!@brief Kawamoto-Schmidt phase*/
  std::valarray<double> kso_; 
  
  std::valarray<size_t> sle_;  /*!@brief generalized slices for e/o */ 
  std::valarray<size_t> slo_;

  void read_bdry(XML::node);

  int esec(int hs)const {return SiteIndex_EvenOdd::instance()->esec(hs);}
  int osec(int hs)const {return SiteIndex_EvenOdd::instance()->osec(hs);}

  void set_ksphase();
  void set_ustag();
  
  void multPoe(Field&,const Field&,int)const;
  void multPeo(Field&,const Field&,int)const;

  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}

public:
  Dirac_staggered_EvenOdd(double mass,const Field* u,
			  const BoundaryCond* bdry = NULL)
    :mq_(mass),u_(u),
     Nvh_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(bdry),
     ff_(Nvh_),gf_(Nvh_),gff_(2*Nvh_),
     ue_(Nvh_),uo_(Nvh_),
     fsize_(ff_.size()),gsize_(gff_.size()),
     kse_(1.0,Nvh_*Ndim_),kso_(1.0,Nvh_*Ndim_),
     sle_(ff_.get_sub(SiteIndex_EvenOdd::instance()->esec())),
     slo_(ff_.get_sub(SiteIndex_EvenOdd::instance()->osec())){
    Mapping::init_shiftField_EvenOdd();
    //
    if(bdry_== NULL) bdry_= new BoundaryCond_periodic;
    set_ksphase();
    set_ustag();
  }

  Dirac_staggered_EvenOdd(const XML::node& stg_node,const Field* u)
   :u_(u),
    Nvh_(CommonPrms::instance()->Nvol()/2),
    Ndim_(CommonPrms::instance()->Ndim()),
    bdry_(NULL),
    ff_(Nvh_),gf_(Nvh_),gff_(2*Nvh_),
    ue_(Nvh_),uo_(Nvh_),
    fsize_(ff_.size()),gsize_(gff_.size()),
    kse_(1.0,Nvh_*Ndim_),kso_(1.0,Nvh_*Ndim_),
    sle_(ff_.get_sub(SiteIndex_EvenOdd::instance()->esec())),
    slo_(ff_.get_sub(SiteIndex_EvenOdd::instance()->osec())){
    Mapping::init_shiftField_EvenOdd();
    //
    XML::read(stg_node,"mass",mq_);
    XML::node bdnode = stg_node;
    XML::descend(bdnode,"boundary_condition",MANDATORY);
    CCIO::cout<<"creating BC"<<std::endl;
    bdry_= createBC(bdnode);
    CCIO::cout<<"BC created"<<std::endl;
    set_ksphase();
    CCIO::cout<<"ksphase created"<<std::endl;
    set_ustag();
    CCIO::cout<<"ustag created"<<std::endl;
  }

  ~Dirac_staggered_EvenOdd(){ if(bdry_!= NULL) delete bdry_;}

  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  
  const Field mult_full(const Field&)const;
  const Field mult_full_dag(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // EvenOdd operator has no preconditioner now 
  const Field mult_prec     (const Field&f)const{return f;}
  const Field mult_dag_prec (const Field&f)const{return f;}
  const Field left_prec     (const Field&f)const{return f;}
  const Field right_prec    (const Field&f)const{return f;}
  const Field left_dag_prec (const Field&f)const{return f;}
  const Field right_dag_prec(const Field&f)const{return f;}
  //////////////////////////////////////////////////////////////
  
  const Field mult_eo(const Field& f) const; 
  const Field mult_oe(const Field& f) const; 
  const Field mult_eo_dag(const Field& f) const;
  const Field mult_oe_dag(const Field& f) const;
  const Field mult_oo(const Field& f)const {return f;}
  const Field mult_ee(const Field& f)const {return f;}
  const Field mult_oo_inv(const Field& f)const {return f;}
  const Field mult_ee_inv(const Field& f)const {return f;}

  const Field md_force(const Field&,const Field&) const;

  void update_internal_state(){ set_ustag(); }
  double get_mq()const{ return mq_;}
  const sfmt_t get_fermionFormat()const{return ff_;}
  const std::vector<int> get_gsite() const{
    return SiteIndex_EvenOdd::instance()->get_gsite(); }
};

#endif
