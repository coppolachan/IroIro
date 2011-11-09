/*!
 * @file qprop_MultiShift.hpp
 *
 * @brief Quark Propagator Qprop_MultiShift class definition
 *
 */

#ifndef QPROP_MULTISHIFT_INCLUDED
#define QPROP_MULTISHIFT_INCLUDED

#include "include/commonPrms.h"
#include "include/pugi_interface.h"
#include "Measurements/FermionicM/source.h"
#include "Solver/multiShiftSolver.h"
#include "Dirac_ops/dirac.h"
#include <vector>

typedef std::vector<Field> prop_t;

class Qprop_MultiShift{
private:
  const Dirac* D_;
  const MultiShiftSolver* slv_;
  int Nc_;
  int Nd_;
  int fsize_; 

public:  
  Qprop_MultiShift(const Dirac* D,
		   const MultiShiftSolver* Solver)
    :D_(D),
     slv_(Solver),
     Nc_(   CommonPrms::instance()->Nc()),
     Nd_(   CommonPrms::instance()->Nd()),
     fsize_(CommonPrms::instance()->Nvol()*Nc_*Nd_*2){}

  ~Qprop_MultiShift(){}

  void calc(prop_t& xq,
	    Source& source,
	    std::vector<double> mass_shifts) const;
};


#endif
