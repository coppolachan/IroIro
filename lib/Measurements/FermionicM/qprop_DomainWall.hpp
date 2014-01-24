/*!
 * @file qprop_DomainWall.hpp
 * @brief Definition of Qprop domain wall fermions QpropDWF
 Time-stamp: <2013-12-05 09:39:59 noaki>
 */
#ifndef QPROP_DOMAINWALL_INCLUDED
#define QPROP_DOMAINWALL_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/dirac_WilsonLike.hpp"
#include "include/commonPrms.hpp"
#include "source.hpp"

typedef std::vector<Field> prop_t;

class QpropDWF : public QuarkPropagator{
  const Dirac_DomainWall_4D& Dgw_;
  int Nc_;
  int Nd_;

public:
  QpropDWF(const Dirac_DomainWall_4D& Ddw4D)
    :Dgw_(Ddw4D),
     Nc_(CommonPrms::instance()->Nc()),
     Nd_(CommonPrms::instance()->Nd()){}

  ~QpropDWF(){}
  
  const Dirac_DomainWall_4D* getKernel() const {return &Dgw_;}
  void calc(prop_t& xq,Source& src) const;
  void calc(prop_t& xq,Source& src, int, int) const;/*!< For tests */
  void calc(prop_t& xq,const prop_t& prp) const;
};

#endif
