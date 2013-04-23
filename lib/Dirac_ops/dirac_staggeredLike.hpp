/*! @file dirac_staggeredLike.hpp
 *  @brief defines the abstract base classes for DiracStaggerdLike 
 * Time-stamp: <2013-04-20 07:23:02 noaki>
 */
#ifndef DIRAC_STAGGEREDLIKE_INCLUDED
#define DIRAC_STAGGEREDLIKE_INCLUDED

#include "dirac.hpp"
#include "include/format_G.h"

typedef Format::Format_G gfmt_t;     // link variables (fundamental rep.)

namespace Dstagg{
  enum Dtype{DdagDee=0,DdagDoo=1,Dfull=2};

  void set_ksphase(std::valarray<double>& ev,std::valarray<double>& od,int Nv);
}

/*!
 * @class DiracStaggeredLike
 * @brief Declaration of abstract base class for Dirac operators of Staggered type
 */
class DiracStaggeredLike : public Dirac {
 public:
  virtual ~DiracStaggeredLike(){}
};

/*!
 * @class DiracStaggeredLike_EvenOdd
 * @brief Declaration of abstract base class for Dirac operators of Staggered type
 */
class DiracStaggeredEvenOddLike : public DiracStaggeredLike {
public:
  virtual ~DiracStaggeredEvenOddLike(){}

  virtual const Field mult_eo(const Field&) const =0;
  virtual const Field mult_oe(const Field&) const =0;
  virtual const Field mult_eo_dag(const Field&) const =0;
  virtual const Field mult_oe_dag(const Field&) const =0;

  virtual const Field mult_ee(const Field&) const =0;
  virtual const Field mult_oo(const Field&) const =0;
  virtual const Field mult_ee_inv(const Field&) const =0;
  virtual const Field mult_oo_inv(const Field&) const =0;
  virtual double getMass() const =0;
};
#endif
