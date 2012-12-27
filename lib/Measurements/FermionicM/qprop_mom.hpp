/*! @file qprop_mom.hpp
  @brief declaration of QpropMom class */
#ifndef QPROP_MOM_INCLUDED
#define QPROP_MOM_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/dirac.hpp"
#include "include/format_F.h"
#include "source.hpp"
#include <vector>

class QpropMom{
private:
  const prop_t& Sq_;
  int Lx_,Ly_,Lz_,Lt_;
  std::vector<std::vector<int> > platt_;
  Format::Format_F fmt_;
  void init_mom();
  void fourier_tr(std::vector<double>& Sp_re,
		  std::vector<double>& Sp_im,
		  const std::vector<int>& p) const;
public:
  QpropMom(const prop_t& Sq)
 :Sq_(Sq),
  Lx_(CommonPrms::instance()->Lx()),
  Ly_(CommonPrms::instance()->Ly()),
  Lz_(CommonPrms::instance()->Lz()),
  Lt_(CommonPrms::instance()->Lt()),
  fmt_(CommonPrms::instance()->Nvol()){ init_mom();}
  
  void output()const;
};

#endif
