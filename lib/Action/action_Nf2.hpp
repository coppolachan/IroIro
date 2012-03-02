/*! 
 * @file action_Nf2.h
 *
 * @brief Declaration of Action_Nf2 class
 */
#ifndef ACTION_NF2_INCLUDED
#define ACTION_NF2_INCLUDED

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/SmartConf.hpp"

/*!
 * @class Action_Nf2
 *
 * @brief Class to calculate TwoFlavors HMC action term
 */
class Action_Nf2 :public Action{
private:
  GaugeField* const u_;      /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;        /*!< Linear solver operator */
  FermionField phi_;
  int fermion_size_;
  bool smeared_;
  SmartConf* SmartField_;
  
  FermionField DdagD_inv(const FermionField& src);
  void attach_smearing(SmartConf*);
public:
  /*!
   * @brief Standard constructor 
   * CG solver is assumed
   */
  Action_Nf2(GaugeField* const GField,
	     DiracWilsonLike* const D, 
	     const Solver* Solv,
	     bool smeared = false,
	     SmartConf* SmearObj = NULL)
    :u_(GField),
     D_(D),
     slv_(Solv),
     smeared_(smeared),
     fermion_size_(D->fsize()){
    if (smeared_ && SmearObj !=NULL) attach_smearing(SmearObj);
  }

  ~Action_Nf2(){}

  void init(const RandNum& rand);
  void observer_update();

  double calc_H();
  GaugeField md_force();

};
#endif
