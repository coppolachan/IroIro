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
  virtual const Field md_force(const Field& eta,const Field& zeta)const = 0;
};

/*
 * @class DiracWilsonLike
 *
 * @brief Declaration of abstract base class for Dirac operators of Wilson type
 */
class DiracWilsonLike : public Dirac {
public:
  virtual ~DiracWilsonLike(){};
  virtual const Field gamma5(const Field&) const = 0;

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


