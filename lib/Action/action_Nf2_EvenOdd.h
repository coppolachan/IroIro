/*! 
 * @file action_Nf2_EvenOdd.h
 *
 * @brief Declaration of Action_Nf2_EvenOdd class
 *
 */
#ifndef ACTION_NF2_EVENODD_INCLUDED
#define ACTION_NF2_EVENODD_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.h"
#include "Dirac_ops/dirac.h"
#include "Solver/solver.h"

/*!
 * @class Action_Nf2_EvenOdd
 *
 * @brief Class to calculate TwoFlavors HMC action term
 *
 */
class Action_Nf2_EvenOdd :public Action{
private:
  Field* const u_; /*!< The gauge field */
  const Dirac* D_; /*!< Dirac Kernel operator */ 
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
  Action_Nf2_EvenOdd(Field* const GField,
		     const Dirac* D, 
		     const Solver* Solv)
    :u_(GField),D_(D),slv_(Solv),
     fsize_(D_->fsize()),phi_(fsize_){}
  
  ~Action_Nf2(){}

  void init(const RandNum& rand,const void* = 0);
  Field md_force(const void* = 0);
  double calc_H();
};
#endif
