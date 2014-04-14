/*!
 * @file qprop_DomainWall.hpp
 * @brief Definition of Qprop domain wall fermions QpropDWF
 Time-stamp: <2014-04-11 14:11:57 noaki>
 */
#ifndef QPROP_DOMAINWALL_INCLUDED
#define QPROP_DOMAINWALL_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/dirac_WilsonLike.hpp"
#include "include/commonPrms.hpp"
#include "source.hpp"

#define PHYSICAL_QUARK_PROPAGATOR

typedef std::vector<Field> prop_t;

class QpropDWF : public QuarkPropagator{
  const Dirac_DomainWall_4D& Dgw_;
  int Nc_;
  int Nd_;

public:
  QpropDWF(const Dirac_DomainWall_4D& Ddw4D)
    :Dgw_(Ddw4D),
     Nc_(CommonPrms::instance()->Nc()),
     Nd_(CommonPrms::instance()->Nd()){
#ifdef PHYSICAL_QUARK_PROPAGATOR
    CCIO::cout<<"QpropDWF: obtaining physical quark propagator\n";
#else 
    CCIO::cout<<"QpropDWF: obtaining 1/Ddwf(4)\n";
#endif
  }

  ~QpropDWF(){}
  
  const Dirac_DomainWall_4D* getKernel() const {return &Dgw_;}
  void calc(prop_t& xq,Source& src) const;
  void calc(prop_t& xq,Source& src, int, int) const;/*!< For tests */
  void calc(prop_t& xq,const prop_t& prp) const;
};

#endif
