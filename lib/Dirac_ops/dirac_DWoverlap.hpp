/*! @file dirac_DWoverlap.hpp
 *  @brief A wrapper of Dirac_DomainWall_4D which acts as the overlap Dirac op.
 * Time-stamp: <2013-12-05 10:00:32 noaki>
 */
#ifndef DIRAC_DWOVERLAP_INCLUDED
#define DIRAC_DWOVERLAP_INCLUDED

#include "dirac_WilsonLike.hpp"
#include "EigenModes/eigenModes.hpp"

class Dirac_DWoverlap :public DiracWilsonLike{

  const Dirac_DomainWall_4D* Ddw4d_;
  const EigenModes* ems_;
  std::vector<double> sgnEv_;
  double mq_;
public:  
  Dirac_DWoverlap(const Dirac_DomainWall_4D* Ddw4d,
		  EigenModes* const ems)
    :Ddw4d_(Ddw4d),ems_(ems),mq_(Ddw4d->getMass()){

    int esize = ems->evals_.size();
    for(int e=0; e<esize; ++e)
      sgnEv_.push_back(ems->evals_[e]/fabs(ems->evals_[e]));
  }

  size_t fsize()const{return Ddw4d_->fsize();}
  size_t gsize()const{return Ddw4d_->gsize();}
  double getMass() const {return mq_;}
  
  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  const Field gamma5  (const Field& f)const{ return Ddw4d_->gamma5(f);}
  const Field* getGaugeField_ptr()const{ return Ddw4d_->getGaugeField_ptr(); }
};
 
#endif
