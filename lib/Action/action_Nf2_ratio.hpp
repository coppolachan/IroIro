/*!
 * @file action_Nf2_ratio.h
 *
 * @brief Declaration of Action_Nf2_ratio class
 *
 */
#ifndef ACTION_RATIO_INCLUDED
#define ACTION_RATIO_INCLUDED

#include "Tools/randNum_MP.h"
#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"


/*!
 * @class Action_Nf2_ratio
 *
 * @brief Class to calculate HMC action term \f$det(D1)/det(D2)\f$
 *
 */
class Action_Nf2_ratio : public Action{
private:
  GaugeField* const u_;
  const DiracWilsonLike* D1_;
  const DiracWilsonLike* D2_;
  const Solver* slv1_;
  const Solver* slv2_;
  FermionField phi_;
  
  FermionField DdagD1_inv(const FermionField& src);
  FermionField DdagD2_inv(const FermionField& src);
  
 public:
  Action_Nf2_ratio(GaugeField* const GField, 
		   const DiracWilsonLike* D1,const DiracWilsonLike* D2,
		   const Solver* Solv1,const Solver* Solv2)
    :u_(GField),
     D1_(D1), D2_(D2),
     slv1_(Solv1), slv2_(Solv2){}
  
  ~Action_Nf2_ratio();
  
  void init(const RandNum& rand);  
  void observer_update(){};

  double calc_H();
  GaugeField md_force();


};

#endif
