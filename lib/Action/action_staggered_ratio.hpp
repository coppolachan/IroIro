/*!
 * @file action_staggered_ratio.h
 * @brief Declaration of Action_staggered_ratio class
 */
#ifndef ACTION_STAGGERED_RATIO_INCLUDED
#define ACTION_STAGGERED_RATIO_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/smartConf.hpp"
/*!
 * @class Action_staggered_ratio
 * @brief Class to calculate HMC action term \f$det(D1)/det(D2)\f$
 */
class Action_staggered_ratio : public Action{
private:
  GaugeField* const u_;
  DiracStaggeredEvenOddLike* const D_;
  const Solver* slv1e_;
  const Solver* slv1o_;
  const Solver* slv2e_;
  const size_t fsize_;
  Field phi_;
  bool smeared_;
  double mr_;
  SmartConf* smart_conf_;

  Field DdagD1e_inv(const Field& src);
  Field DdagD1o_inv(const Field& src);
  Field DdagD2e_inv(const Field& src);
  void attach_smearing(SmartConf*);
 public:
  Action_staggered_ratio(GaugeField* const GField, 
			 DiracStaggeredEvenOddLike* const D1,
			 DiracStaggeredEvenOddLike* const D2,
			 const Solver* Solv1e,const Solver* Solv1o,
			 const Solver* Solv2e,
			 bool smeared = false,
			 SmartConf* smart_conf = NULL)
    :u_(GField),D_(D1),
     slv1e_(Solv1e),slv1o_(Solv1o),slv2e_(Solv2e),
     fsize_(D1->fsize()),
     phi_(fsize_),mr_(D1->get_mq()/D2->get_mq()),
     smeared_(smeared){
    if(smeared_ && smart_conf !=NULL) attach_smearing(smart_conf);
  }
  
  ~Action_staggered_ratio(){}
  
  void init(const RandNum& rand);  
  double calc_H();
  GaugeField md_force();
  void observer_update();

};

#endif
