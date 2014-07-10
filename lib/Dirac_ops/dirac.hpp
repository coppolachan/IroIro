/*!
 * @file dirac.hpp
 * @brief Defines the abstract base classs for Dirac operators
 * Time-stamp: <2014-07-08 16:40:39 noaki>
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
  virtual double getMass()const = 0;

  virtual const Field* getGaugeField_ptr()const =0;
  virtual const Field mult    (const Field&)const = 0;
  virtual const Field mult_dag(const Field&)const = 0;

  virtual const Field md_force(const Field& eta,const Field& zeta)const{}
  virtual void update_internal_state(){}
};

#endif


