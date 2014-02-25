/*! 
  @file chiral_condensate.hpp
  @brief Declaration of Chiral condensate measurement class ChiralCond
 */
#ifndef CHIRAL_COND_INCLUDED
#define CHIRAL_COND_INCLUDED
 
#include "chiral_condensate_abs.hpp"
#include "include/errors.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Measurements/FermionicM/source.hpp"
#include "Solver/solver_Factory.hpp"

/*!
 * @class ChiralCondStd
 * @brief Calculates the Chiral condensate \f$ \rangle \bar\psi \psi \f$
 * 
 * It is designed for any operator except the DWF
 *
 */
class ChiralCondStd : public ChiralCondensate{
private:
  const Dirac* D_;/*!< @brief %Dirac operator */
  const Solver* slv_;/*!< @brief %Solver for the inversion */

  Field invert(Field&)const; 
  int fsize()const;
public:
  /*! @brief Public constructor */
  ChiralCondStd(const Dirac* D,
		const Solver* Solver)
    :D_(D),
     slv_(Solver){
    if (!(slv_->check_DdagD())) {
      ErrorString msg;
      msg << "Input Solver has no Fopr_DdagD Kernel.\n";
      Errors::BaseErr("General Error", msg);
    }
  }
  /*! @brief Public destructor */
  ~ChiralCondStd(){}
};

/*!
 * @class ChiralCondDWF
 * @brief Calculates the Chiral condensate \f$ \rangle \bar\psi \psi \f$
 * 
 * It is designed for the DWF
 *
 */
class ChiralCondDWF : public ChiralCondensate{
private:
  const Dirac_DomainWall_4D& Ddw_;/*!< @brief %Dirac operator */
  Field invert(Field&)const; 
  int fsize()const;
  double one_minus_m_inv; // (1-m)^-1
public:
  ChiralCondDWF(const Dirac_DomainWall_4D& Ddw4D)
    :Ddw_(Ddw4D),one_minus_m_inv(1.0/(1.0-Ddw_.getMass())){}

  /*! @brief Public destructor */
  ~ChiralCondDWF(){}
};


#endif
