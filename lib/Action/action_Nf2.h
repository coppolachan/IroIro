/*! 
 * @file action_Nf2.h
 *
 * @brief Declaration of Action_Nf2 class
 *
 */
#ifndef ACTION_PSF_INCLUDED
#define ACTION_PSF_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.h"
#include "include/fopr.h"
#include "include/format_G.h"
#include "Dirac_ops/dirac_Operator_Factory.hpp"
#include "Solver/solver_Factory.hpp"


/*!
 * @class Action_Nf2
 *
 * @brief Class to calculate TwoFlavors HMC action term
 *
 */
class Action_Nf2 :public Action{
private:
  Field* const u_; /*!< The gauge field */
  const Dirac* D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;  /*!< Linear solver operator */

  Field phi_;
  size_t fsize_;
  
  Field DdagD_inv(const Field& src);
  
public:
  /*!
   * @brief Standard constructor 
   * 
   * CG solver is assumed
   */
  Action_Nf2(Field* const GField,
	     const Dirac* D, 
	     const Solver* Solv)
    :u_(    GField),
     D_(D),
     slv_(Solv),
     fsize_(D_->fsize()){}

  ~Action_Nf2(){}
  

  void init(Field&,const RandNum& rand,const void* = 0);
  
  void init(Field& P,const RandNum& rand,const Field& U, const void* = 0);

  double calc_H();

  Field md_force(const void* = 0);

  Field md_force(const Field& U, const void* = 0);
};
#endif
