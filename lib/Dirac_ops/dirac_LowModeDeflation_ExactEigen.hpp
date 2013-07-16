/*!@file  dirac_LowModeDeflation_ExactEigen.hpp
 * @brief concrete class for low-mode preconditioning with exact low-modes
 *Time-stamp: <2013-07-02 14:37:05 noaki>
 */
#ifndef DIRAC_LOWMODEDEFLATION_EXACTEIGEN_INCLUDED
#define DIRAC_LOWMODEDEFLATION_EXACTEIGEN_INCLUDED

#include "field.h"
#include "dirac_LowModeDeflation.hpp"
#include "EigenModes/eigenModes.hpp"

class Dirac_LowModeDeflation_ExactEigen : public Dirac_LowModeDeflation{

  const DiracWilsonLike* Dw_;
  std::vector<Field>* const evecs_;
  std::vector<double>* const evals_;
  std::vector<double> ievals_;
  std::vector<double> evals_conj_;
  std::vector<double> ievals_conj_;

  void setEigenArrays();
public:
  Dirac_LowModeDeflation_ExactEigen(const DiracWilsonLike* Dw,
				    EigenModes* const ems)
    :Dw_(Dw),evecs_(&(ems->evecs_)),evals_(&(ems->evals_)){setEigenArrays();}

  size_t fsize()const {return Dw_->fsize();}
  size_t gsize()const {return Dw_->gsize();}
  double getMass()const{return Dw_->getMass();}  

  const Field* getGaugeField_ptr()const{ return Dw_->getGaugeField_ptr(); }

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field mult_little(const Field&)const;
  const Field mult_little_inv(const Field&)const;
  const Field mult_little_dag(const Field&)const;
  const Field mult_little_inv_dag(const Field&)const;
  
  const Field mult_left(const Field&)const;
  const Field mult_right(const Field& f)const{return f;}

  const Field gamma5(const Field& f)const{ return Dw_->gamma5(f);}
};

#endif
