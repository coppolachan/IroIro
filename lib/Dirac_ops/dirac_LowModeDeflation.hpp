/*!@file  dirac_LowModeDeflation.hpp
 * @brief abstract class for low-mode preconditioning
 * Time-stamp: <2013-07-02 10:42:22 noaki>
 */
#ifndef DIRAC_LOWMODEDEFLATION_INCLUDED
#define DIRAC_LOWMODEDEFLATION_INCLUDED

#include "dirac_WilsonLike.hpp"

class Dirac_LowModeDeflation : public DiracWilsonLike{
public:
  virtual const Field mult_little(const Field&)const = 0;
  virtual const Field mult_little_inv(const Field&)const = 0;
  virtual const Field mult_little_dag(const Field&)const = 0;
  virtual const Field mult_little_inv_dag(const Field&)const = 0;
  
  virtual const Field mult_left(const Field&)const = 0;
  virtual const Field mult_right(const Field&)const = 0;
  
  virtual ~Dirac_LowModeDeflation(){}
};

#endif
