/*! 
 * @file action_Nf2.h
 *
 * @brief Declaration of Action_Nf2 class
 *
 */
#ifndef ACTION_NF2_INCLUDED
#define ACTION_NF2_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.hpp"
#include "Dirac_ops/dirac.h"
#include "Solver/solver.hpp"

/*!
 * @class Action_Nf2
 *
 * @brief Class to calculate TwoFlavors HMC action term
 *
 */
class Action_Nf2 :public Action{
private:
  Field* const u_; /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;  /*!< Linear solver operator */
  size_t fsize_;
  Field phi_;
  
  Field DdagD_inv(const Field& src);
  
public:
  /*!
   * @brief Standard constructor 
   * 
   * CG solver is assumed
   */
  Action_Nf2(Field* const GField,
	     DiracWilsonLike* const D, 
	     const Solver* Solv)
    :u_(GField),D_(D),slv_(Solv),
     fsize_(D_->fsize()),phi_(fsize_){}

  ~Action_Nf2(){}

  void init(const RandNum& rand,const void* = 0);
  Field md_force(const void* = 0);
  double calc_H();
  void observer_update();
};
#endif
