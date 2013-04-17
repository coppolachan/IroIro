/*! 
 * @file action_staggered.hpp
 * @brief Declaration of Action_staggered class
 */
#ifndef ACTION_STAGGERED_INCLUDED
#define ACTION_STAGGERED_INCLUDED

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/smartConf.hpp"

class Action_staggered :public Action{
private:
  GaugeField* const u_;                 /*!< The gauge field */
  DiracStaggeredEvenOddLike* const D_;  /*!< Dirac Kernel operator */ 
  const Solver* slv_;                   /*!< Linear solver operator */
  const size_t fsize_;  
  Field phi_;
  bool smeared_;
  SmartConf* smart_conf_;
  
  Field DdagD_inv(const Field& src);
  void attach_smearing(SmartConf*);
 
public:
  /*! @brief Standard constructor */
  Action_staggered(GaugeField* const GField,
		   DiracStaggeredEvenOddLike* const D, 
		   const Solver* Solv,
		   bool smeared = false,
		   SmartConf* smart_conf = NULL)
    :u_(GField),D_(D),slv_(Solv),
     fsize_(D->fsize()),phi_(fsize_),
     smeared_(smeared){
     if (smeared_ && smart_conf !=NULL) attach_smearing(smart_conf);
  }

  ~Action_staggered(){}

  void init(const RandNum& rand);
  double calc_H();
  GaugeField md_force();
  void observer_update();

};
#endif

