/*! 
  @file chiral_condensate.hpp
  @brief Declaration of Chiral condensate measurement class ChiralCond
 */
#ifndef CHIRAL_COND_INCLUDED
#define CHIRAL_COND_INCLUDED
 
#include "chiral_condensate_abs.hpp"
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

public:
  /*! @brief Public constructor */
  ChiralCondStd(const Dirac* D,
	const Solver* Solver)
    :D_(D),
     slv_(Solver){
    if (!(slv_->check_DdagD())) {
      std::cerr<< "Input Solver has no Fopr_DdagD Kernel. ";
      abort();
    }
  }
  /*! @brief Public destructor */
  ~ChiralCondStd(){}
  
  /*! @brief Calculates the chiral condensate
   */
  double calc(Source&) const;
};

/*!
 * @class ChiralCondStd
 * @brief Calculates the Chiral condensate \f$ \rangle \bar\psi \psi \f$
 * 
 * It is designed for the DWF
 *
 */
class ChiralCondDWF : public ChiralCondensate{
private:
  const Dirac_DomainWall_4D& Dgw_;/*!< @brief %Dirac operator */
  int Nc_;
  int Nd_;

public:
  ChiralCondDWF(const Dirac_DomainWall_4D& Ddw4D)
    :Dgw_(Ddw4D),
     Nc_(CommonPrms::instance()->Nc()),
     Nd_(CommonPrms::instance()->Nd()){}

  /*! @brief Public destructor */
  ~ChiralCondDWF(){}
  
  /*! @brief Calculates the chiral condensate
   */
  double calc(Source&) const;
};



#endif
