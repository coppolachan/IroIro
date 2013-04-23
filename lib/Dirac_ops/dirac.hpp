/*!
 * @file dirac.hpp
 * @brief Defines the abstract base class for Dirac operators
 */
#ifndef DIRAC_H_
#define DIRAC_H_

#include "include/field.h"
#include "include/pugi_interface.h"

enum {Op, Dag};

namespace Dop{
  // for Wilson_EvenOdd
  struct EOtag{};    
  struct OEtag{};    

  double read_mass(const XML::node& node);
}
class RandNum;

/*
 *! @class Dirac
 * @brief Declaration of abstract base class for Dirac operators
 */
class Dirac {
public:
  virtual ~Dirac(){}
  virtual size_t fsize() const = 0;
  virtual size_t gsize() const = 0;
  virtual int Nvol()const = 0;

  virtual const Field mult    (const Field&)const = 0;
  virtual const Field mult_dag(const Field&)const = 0;

  /////////////////////////////////Preconditioned versions
  virtual const Field mult_prec    (const Field&)const = 0;
  virtual const Field mult_dag_prec(const Field&)const = 0;

  virtual const Field left_prec     (const Field&)const = 0;
  virtual const Field right_prec    (const Field&)const = 0;
  virtual const Field left_dag_prec (const Field&)const = 0;
  virtual const Field right_dag_prec(const Field&)const = 0;
  ////////////////////////////////////////////

  virtual const Field md_force(const Field& eta,const Field& zeta)const = 0;
  virtual void update_internal_state() = 0;
  virtual void get_RandGauss(std::valarray<double>&,const RandNum&)const =0;
};

/*!
 * @class DiracWilsonLike
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike : public Dirac {
public:
  virtual ~DiracWilsonLike(){}
  virtual const Field gamma5(const Field&) const = 0;
};

/*!
 * @class DiracWilsonLike_EvenOdd
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike_EvenOdd : public DiracWilsonLike {
public:
  virtual ~DiracWilsonLike_EvenOdd(){}

  virtual const Field mult_eo(const Field&) const =0;
  virtual const Field mult_oe(const Field&) const =0;
  virtual const Field mult_eo_dag(const Field&) const =0;
  virtual const Field mult_oe_dag(const Field&) const =0;

  virtual const Field mult_ee(const Field&) const =0;
  virtual const Field mult_oo(const Field&) const =0;
  virtual const Field mult_ee_inv(const Field&) const =0;
  virtual const Field mult_oo_inv(const Field&) const =0;
};

/*!
 * @class Dirac_OptimalDomainWall_4D
 * @brief Declaration of abstract base class for 4D-reducted Domain-Wall fermions
 */
class Dirac_optimalDomainWall_4D : public DiracWilsonLike {
public:
  virtual ~Dirac_optimalDomainWall_4D(){}

  virtual double getMass()const =0;
  virtual const Field mult_inv(const Field&) const =0;
  virtual const Field mult_dag_inv(const Field&) const =0;
};

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


