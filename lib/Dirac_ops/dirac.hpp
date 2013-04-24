/*!
 * @file dirac.hpp
 * @brief Defines the abstract base classs for Dirac operators
 * Time-stamp: <2013-04-24 09:54:26 noaki>
 */
#ifndef DIRAC_H_
#define DIRAC_H_

#include "include/field.h"
#include "include/pugi_interface.h"

/*
 *! @class Dirac
 * @brief Declaration of abstract base class for Dirac operators
 */
class Dirac {
public:
  virtual ~Dirac(){}
  virtual size_t fsize() const = 0;
  virtual size_t gsize() const = 0;
  //  virtual int Nvol()const = 0;

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

  virtual const Field md_force(const Field& eta,const Field& zeta)const{}
  virtual void update_internal_state() = 0;
};

#endif


