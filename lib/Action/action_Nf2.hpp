/*! 
 * @file action_Nf2.hpp
 * @brief Declaration of the Action_Nf2 class
 *
 * Time-stamp: <2013-04-22 17:04:47 neo>
 */
#ifndef ACTION_NF2_INCLUDED
#define ACTION_NF2_INCLUDED

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/smartConf.hpp"

#include "Dirac_ops/boundaryCond.hpp"

/*!
 * @class Action_Nf
 * @brief Class to calculate 2 Flavors of Dirac fermions
 */
class Action_Nf2 :public Action{
private:
  GaugeField* const u_;      /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;        /*!< Linear solver operator */
  const size_t fsize_;       /*!< Fermion field size */
  Field phi_;                /*!< Pseudofermion field */
  bool smeared_;             /*!< Asserts if the action is smeared */
  SmartConf* smart_conf_;    /*!< Pointer to the Smart Configurations */
  
  BoundaryCond* BC;

  Field DdagD_inv(const Field& src);
  void attach_smearing(SmartConf*);
public:
  /*! @brief Standard constructor */
  Action_Nf2(GaugeField* const GField,
	     DiracWilsonLike* const D, 
	     const Solver* Solv,
	     bool smeared = false,
	     SmartConf* smart_conf = NULL)
    :u_(GField),
     D_(D),
     slv_(Solv),
     fsize_(D->fsize()),
     phi_(fsize_),
     smeared_(smeared),
     BC(new BoundaryCond_antiPeriodic(TDIR)){
     if (smeared_ && smart_conf !=NULL) attach_smearing(smart_conf);
  }

  ~Action_Nf2(){}

 /*! @brief Initializes the pseudofermion fields */ 
  void init(const RandNum& rand);

  /*! @brief Calculates the action contribution to H */
  double calc_H();

  /*! @brief Force term for the HMC update */
  GaugeField md_force();

  void observer_update();
};
#endif
