/*!
 * @file qprop_optimalDomainWall.hpp
 *
 * @brief Definition of Qprop optimal domain wall fermions QpropDWF
 *
 */
#ifndef QPROP_OPTIMALDOMAINWALL_INCLUDED
#define QPROP_OPTIMALDOMAINWALL_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/dirac_DomainWall_4D.hpp"
#include "source.h"

typedef std::vector<Field> prop_t;

class QpropDWF : public QuarkPropagator
{
  const Dirac_optimalDomainWall_4D Dgw_;
  int Nc_;
  int Nd_;

public:
  QpropDWF(Dirac_optimalDomainWall_4D& DWF_Kernel)
    :Dgw_(DWF_Kernel),
     Nc_(   CommonPrms::instance()->Nc()),
     Nd_(   CommonPrms::instance()->Nd()){}
  
  QpropDWF(const Dirac_optimalDomainWall& DWF_5dKernel,
	   double scnd = 1.0e-14,
	   int Niter =500)
    :Dgw_(DWF_5dKernel,scnd,scnd,Niter),
     Nc_(   CommonPrms::instance()->Nc()),
     Nd_(   CommonPrms::instance()->Nd()){}

  ~QpropDWF(){}
  
  void calc(prop_t& xq,Source& src) const;
  void calc(prop_t& xq,Source& src, int, int) const;/*!< For tests */
  void calc(prop_t& xq,const prop_t& prp) const;
};

#endif
