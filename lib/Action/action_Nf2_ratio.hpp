/*!
 * @file action_Nf2_ratio.hpp
 * @brief Declaration of Action_Nf2_ratio class
 Time-stamp: <2013-07-16 17:25:40 cossu>
 */
#ifndef ACTION_NF2_RATIO_INCLUDED
#define ACTION_NF2_RATIO_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.hpp"
#include "Dirac_ops/dirac_WilsonLike.hpp"
#include "Solver/solver.hpp"
#include "Smearing/smartConf.hpp"
/*!
 * @class Action_Nf2_ratio
 * @brief Class to calculate HMC action term \f$det(D1)/det(D2)\f$
 */
class Action_Nf2_ratio : public Action{
private:
  GaugeField* const u_;
  DiracWilsonLike* D1_;
  DiracWilsonLike* D2_;
  const Solver* slv1_;
  const Solver* slv2_;
  const size_t fsize_;
  Field phi_;
  bool smeared_;
  const char* name_;
  SmartConf* smart_conf_;

  Field DdagD1_inv(const Field& src);
  Field DdagD2_inv(const Field& src);
  void attach_smearing(SmartConf*);
 public:
  Action_Nf2_ratio(GaugeField* const GField, 
		   DiracWilsonLike* const D1,DiracWilsonLike* const D2,
		   const Solver* Solv1,const Solver* Solv2,
		   const char* n = "Action_Nf2_ratio",
		   bool smeared = false,
		   SmartConf* smart_conf = NULL)
    :u_(GField),
     D1_(D1), D2_(D2),
     slv1_(Solv1), slv2_(Solv2),
     fsize_(D1->fsize()),
     phi_(fsize_),
     name_(n),
     smeared_(smeared){
    if (smeared_ && smart_conf !=NULL) attach_smearing(smart_conf);
  }
  
  ~Action_Nf2_ratio(){}
  
  void init(const RandNum& rand);  
  void observer_update();
  double calc_H();
  GaugeField md_force();
};

#endif
