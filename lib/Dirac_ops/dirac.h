/*
 * @file dirac.h
 *
 * @brief Defines the abstract base classs for Dirac operators
 *
 */

#ifndef DIRAC_H_
#define DIRAC_H_

#include "include/field.h"

enum {Op, Dag};


/*
 * @class Dirac
 *
 * @brief Declaration of abstract base class for Dirac operators
 */
class Dirac {
public:
  virtual ~Dirac(){};
  virtual size_t fsize() const = 0;
  virtual size_t gsize() const = 0;

  virtual const Field operator()(int OpType, const Field&)const = 0;


  virtual const Field mult(const Field&)const = 0;
  virtual const Field mult_dag(const Field&)const = 0;
  //Preconditioned versions
  virtual const Field mult_prec(const Field&)const = 0;
  virtual const Field mult_dag_prec(const Field&)const = 0;

  virtual const Field md_force(const Field& eta,const Field& zeta)const = 0;
};

/*
 * @class DiracWilsonLike
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike : public Dirac {
public:
  virtual ~DiracWilsonLike(){};
  virtual const Field gamma5(const Field&) const = 0;
};

/*
 * @class DiracWilsonLike_EvenOdd
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike_EvenOdd : public Dirac {
public:
  virtual ~DiracWilsonLike_EvenOdd(){};
  virtual const Field gamma5(const Field&) const = 0;

  virtual const Field mult_eo(const Field&) const =0;
  virtual const Field mult_oe(const Field&) const =0;
  virtual const Field mult_eo_dag(const Field&) const =0;
  virtual const Field mult_oe_dag(const Field&) const =0;

  virtual const Field mult_ee(const Field&) const =0;
  virtual const Field mult_oo(const Field&) const =0;
  virtual const Field mult_ee_inv(const Field&) const =0;
  virtual const Field mult_oo_inv(const Field&) const =0;
};


/*
 * @class DiracStaggeredLike
 *
 * @brief Declaration of abstract base class for Dirac operators of Staggered type
 */
class DiracStaggeredLike : public Dirac {
 public:
    virtual ~DiracStaggeredLike(){};

};




#endif


