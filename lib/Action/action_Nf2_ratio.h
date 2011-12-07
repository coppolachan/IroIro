/*!
 * @file action_Nf2_ratio.h
 *
 * @brief Declaration of Action_Nf2_ratio class
 *
 */
#ifndef ACTION_RATIO_INCLUDED
#define ACTION_RATIO_INCLUDED

#ifndef RANDNUM_MP_INCLUDED
#include "Tools/randNum_MP.h"
#endif

#ifndef ACTION_INCLUDED
#include "Action/action.h"
#endif

#ifndef FOPR_INCLUDED
#include "include/fopr.h"
#endif

#ifndef FORMAT_G_INCLUDED
#include "include/format_G.h"
#endif

#ifndef FORMAT_F_INCLUDED
#include "include/format_F.h"
#endif

/*!
 * @class Action_Nf2_ratio
 *
 * @brief Class to calculate HMC action term \f$det(D1)/det(D2)\f$
 *
 */
class Action_Nf2_ratio : public Action{
private:
  Field* const u_;
  const Dirac* D1_;
  const Dirac* D2_;
  const Solver* slv1_;
  const Solver* slv2_;
  size_t fsize_;
  Field phi_;
  
  Field DdagD1_inv(const Field& src);
  Field DdagD2_inv(const Field& src);
  
 public:
  Action_Nf2_ratio(Field* const GField, 
		   const Dirac* D1,const Dirac* D2,
		   const Solver* Solv1,const Solver* Solv2)
    :u_(GField),
     D1_(D1), D2_(D2),
     slv1_(Solv1), slv2_(Solv2),
     fsize_(D1_->fsize()),
     phi_(fsize_){}
  
  ~Action_Nf2_ratio();
  
  void init(Field&,const RandNum& rand,const void* = 0);  
  Field md_force(const void* = 0);
  double calc_H();
};

#endif
