/*!@file  dirac_LowModeDeflation_Approx.hpp
 * @brief concrete class for low-mode preconditioning with given subspace
 *Time-stamp: <2013-07-03 10:58:34 noaki>
 */
#ifndef DIRAC_LOWMODEDEFLATION_APPROX_INCLUDED
#define DIRAC_LOWMODEDEFLATION_APPROX_INCLUDED

#include "field.h"
#include "dirac_LowModeDeflation.hpp"
#include "EigenModes/eigenModes.hpp"

class Dirac_LowModeDeflation_Approx : public Dirac_LowModeDeflation{

  const DiracWilsonLike* Dw_;
  const vector_Field* subVecs_;
  
  std::valarray<double> Dsub_;
  std::valarray<double> Dsub_inv_;
  std::valarray<double> Dsub_dag_;
  std::valarray<double> Dsub_dag_inv_;
  
  int Nsub_;
  void rightPrec(Field&,const Field&)const;
  void lowModeProj(Field& w,const Field& f,const std::valarray<double>&)const;
  void calcLittleOp();

public:
  Dirac_LowModeDeflation_Approx(const DiracWilsonLike* Dw,
				const vector_Field* vecs)
    :Dw_(Dw),subVecs_(vecs),Nsub_(vecs->size()){}

  Dirac_LowModeDeflation_Approx(const DiracWilsonLike* Dw,
				EigenModes* const ems)
    :Dw_(Dw),subVecs_(&(ems->evecs_)){}

  size_t fsize() const {return Dw_->fsize();}
  size_t gsize()const{return Dw_->gsize();}
  double getMass() const {return Dw_->getMass();}  

  const Field* getGaugeField_ptr()const{ return Dw_->getGaugeField_ptr(); }

  const Field mult    (const Field&)const;
  const Field mult_dag(const Field&)const;

  const Field mult_little(const Field&)const;
  const Field mult_little_inv(const Field&)const;
  const Field mult_little_dag(const Field&)const;
  const Field mult_little_inv_dag(const Field&)const;
  
  const Field mult_left(const Field&)const;
  const Field mult_right(const Field&)const;

  const Field gamma5(const Field& f)const{ return Dw_->gamma5(f);}
};

#endif
