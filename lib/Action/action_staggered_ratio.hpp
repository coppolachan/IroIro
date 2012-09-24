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
  DiracStaggeredEvenOddLike* const D1_;
  DiracStaggeredEvenOddLike* const D2_;
  const Solver* slv1_;
  const Solver* slv2_;
  const size_t fsize_;
  Field phi_;
  bool smeared_;
  SmartConf* smart_conf_;

  Field DdagD1_inv(const Field& src);
  Field DdagD2_inv(const Field& src);
  void attach_smearing(SmartConf*);
 public:
  Action_staggered_ratio(GaugeField* const GField, 
			 DiracStaggeredEvenOddLike* const D1,
			 DiracStaggeredEvenOddLike* const D2,
			 const Solver* Solv1,const Solver* Solv2,
			 bool smeared = false,
			 SmartConf* smart_conf = NULL)
    :u_(GField),
     D1_(D1), D2_(D2),
     slv1_(Solv1), slv2_(Solv2),
     fsize_(D1->fsize()),
     phi_(fsize_),
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
