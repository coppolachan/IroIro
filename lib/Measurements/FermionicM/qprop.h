/*!
 *
 * @file qprop.h
 *
 * @brief Declaration of Quark propagator class Qprop
 *
 */
#ifndef QPROP_INCLUDED
#define QPROP_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/dirac.h"
#include "Solver/solver.h"
#include "Measurements/FermionicM/source.h"
#include "Dirac_ops/dirac_Factory.hpp"
#include "Solver/solver_Factory.hpp"



/*!
 * @class Qprop
 *
 * @brief Calculates the Quark propagator \f$D^{-1}\f$
 *
 * \f$(D^\dagger D)^{-1}\f$ is applyed to the vector 
 * \f$D^\dagger \psi\f$
 * 
 * \f[ \chi = D^{-1} \psi = (D^\dagger D)^{-1} D^\dagger \psi \f]
 *
 */
class Qprop : public QuarkPropagator{
private:
  const Dirac* D_;/*!< @brief %Dirac operator needed to get \f$(D^\dagger)^{-1}\f$ */
  const Solver* slv_;/*!< @brief %Solver for the inversion of \f$(D^\dagger D)\f$ */

  int Nc_;/*!< @brief Number of colors */
  int Nd_;/*!< @brief Number of %Dirac indexes */
  int fsize_;/*!< @brief Fermion field size */ 

public:
  /*! @brief Public constructor */
  Qprop(const Dirac* D,
	const Solver* Solver)
    :D_(D),
     slv_(Solver), 
     Nc_(   CommonPrms::instance()->Nc()),
     Nd_(   CommonPrms::instance()->Nd()),
     fsize_(D_->fsize()){
  
    if (!(slv_->check_DdagD())) {
      std::cerr<< "Input Solver has no Fopr_DdagD Kernel. ";
      abort();
    }
    
  }
  /*! @brief Public destructor */
  ~Qprop(){}
  
  /*! @brief Calculates the propagator 
   *
   * N.B. A typedef std::vector<Field> prop_t;
   * is defined in quark_propagators.hpp
   */
  void calc(prop_t& xq,Source& src) const;
  void calc(prop_t& xq,Source& src, int, int) const;
  void calc(prop_t& xq,const prop_t& prp) const;
};
#endif
