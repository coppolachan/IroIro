//-------------------------------------------------------------
/*!@file   dirac_staggered_EvenOdd.hpp
 * @brief  Definition of the staggered operator
 */
//-------------------------------------------------------------

#ifndef DIRAC_STAGGERED_EVENODD_INCLUDED
#define DIRAC_STAGGERED_EVENODD_INCLUDED

#include "dirac_staggeredLike.hpp"
#include "include/format_S.h"
#include "Main/Geometry/shiftField.hpp"
#include "boundaryCond.hpp"
#include <memory>

typedef Format::Format_S sfmt_t;     // fermions (fundamental rep.)

class Dirac_staggered_EvenOdd :public DiracStaggeredEvenOddLike{
private:
  const Field* const u_;
  int Nvh_,Ndim_;
  double mq_;
  GaugeField ue_,uo_;

  const sfmt_t ff_;  // for all pseudo-fermions 
  const gfmt_t gf_;  // for ue_ and uo_
  const gfmt_t gff_; // for u_ (always full site)

  const size_t fsize_;
  const size_t gsize_;

  const BoundaryCond* bdry_; 
  
  std::valarray<double> kse_;  /*!@brief Kawamoto-Schmidt phase*/
  std::valarray<double> kso_; 
  
  std::valarray<size_t> sle_;  /*!@brief generalized slices for e/o */ 
  std::valarray<size_t> slo_;

  void set_ustag();

  void multPoe(Field&,const Field&,int)const;
  void multPeo(Field&,const Field&,int)const;

  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}
  
  void mult_DdagDee(Field&,const Field&)const;
  void mult_DdagDoo(Field&,const Field&)const;
  void mult_Dfull(    Field&,const Field&)const;
  void mult_Dfull_dag(Field&,const Field&)const;

  static void (Dirac_staggered_EvenOdd::*mult_type[])(Field&,const Field&)const;
  static void (Dirac_staggered_EvenOdd::*mult_dag_type[])(Field&,const Field&)const;
  
  void (Dirac_staggered_EvenOdd::*mult_core)(Field&,const Field&)const;
  void (Dirac_staggered_EvenOdd::*mult_dag_core)(Field&,const Field&)const;

public:
  Dirac_staggered_EvenOdd(double mass,Dstagg::Dtype D,const Field* u,
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
     slo_(ff_.get_sub(SiteIndex_EvenOdd::instance()->osec())),
     mult_core(Dirac_staggered_EvenOdd::mult_type[D]),
     mult_dag_core(Dirac_staggered_EvenOdd::mult_dag_type[D])
  {
    Mapping::init_shiftField_EvenOdd();
    //
    if(bdry_== NULL) bdry_= new BoundaryCond_periodic;
    Dstagg::set_ksphase(kse_,kso_,Nvh_);
    set_ustag();
  }
  
  Dirac_staggered_EvenOdd(const XML::node& stg_node,
			  Dstagg::Dtype D,const Field* u)
    :u_(u),
     Nvh_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     bdry_(NULL),
     ff_(Nvh_),gf_(Nvh_),gff_(2*Nvh_),
     ue_(Nvh_),uo_(Nvh_),
     fsize_(ff_.size()),gsize_(gff_.size()),
     kse_(1.0,Nvh_*Ndim_),kso_(1.0,Nvh_*Ndim_),
     sle_(ff_.get_sub(SiteIndex_EvenOdd::instance()->esec())),
     slo_(ff_.get_sub(SiteIndex_EvenOdd::instance()->osec())),
     mult_core(Dirac_staggered_EvenOdd::mult_type[D]),
     mult_dag_core(Dirac_staggered_EvenOdd::mult_dag_type[D])
  {
    Mapping::init_shiftField_EvenOdd();                                        
    //
    XML::read(stg_node,"mass",mq_);
    XML::node bdnode = stg_node;
    XML::descend(bdnode,"BoundaryCondition",MANDATORY);
    bdry_= createBC(bdnode);
    Dstagg::set_ksphase(kse_,kso_,Nvh_);
    set_ustag();
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
  double getMass()const{ return mq_;}
  int Nvol()const{return ff_.Nvol();}
};

#endif



