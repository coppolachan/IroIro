/*! 
 * @file action_Nf2.h
 * @brief Declaration of Action_Nf2 class
 */
#ifndef ACTION_NF2_INCLUDED
#define ACTION_NF2_INCLUDED

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/smartConf.hpp"

class Action_Nf2 :public Action{
private:
  GaugeField* const u_;      /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;        /*!< Linear solver operator */
  const size_t fsize_;  
  Field phi_;
  SmartConf* smart_conf_;
  
  Field DdagD_inv(const Field& src);
  void attach_smearing(SmartConf*);
public:
  /*! @brief Standard constructor */
  Action_Nf2(GaugeField* const GField,
	     DiracWilsonLike* const D, 
	     const Solver* Solv,
	     SmartConf* smart_conf = NULL)
    :u_(GField),
     D_(D),
     slv_(Solv),
     fsize_(D->fsize()),
     phi_(fsize_),
     smart_conf_(smart_conf){
    if(smart_conf_) 
      assert(u_== smart_conf_->get_current_conf());
  }

  ~Action_Nf2(){}

  void init(const RandNum& rand);
  double calc_H();
  GaugeField md_force();
  void observer_update(){D_->update_internal_state(); }

};
#endif
