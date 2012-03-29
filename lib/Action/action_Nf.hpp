/*! 
 * @file action_Nf2.h
 *
 * @brief Declaration of Action_Nf class
 *
 * Any number of flavours
 */
#ifndef ACTION_NF_INCLUDED
#define ACTION_NF_INCLUDED

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Smearing/SmartConf.hpp"
#include "Tools/RationalApprox/rationalapprox.hpp"


/*!
 * @class Action_Nf
 *
 * @brief Class to calculate N-Flavors RHMC action term
 */
class Action_Nf :public Action{
private:
  GaugeField* const u_;      /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  const Solver* slv_;        /*!< Multishift Linear solver operator */
  FermionField phi_;
  int fermion_size_;
  int n_pseudof_; /*!< Number of pseudofermion fields */
  int n_flav_;    /*!< Number of flavors */
  bool smeared_;
  SmartConf* SmartField_;
  
  // Rational approximations
  RationalApprox PseudoFermionsApprox_;
  RationalApprox MolecularDynApprox_;
  RationalApprox MetropolisApprox_;
 

  FermionField DdagD_inv(const FermionField& src);
  void attach_smearing(SmartConf*);
public:
  /*!
   * @brief Standard constructor 
   * CG solver is assumed
   */
  Action_Nf(GaugeField* const GField,
	    DiracWilsonLike* const D, 
	    const Solver* Solv,
	    int npseu
	    bool smeared = false,
	    SmartConf* SmearObj = NULL)
    :u_(GField),
     D_(D),
     slv_(Solv),
     smeared_(smeared),
     fermion_size_(D->fsize()){
    phi_.resize(fermion_size_); //takes care of EvenOdd and 5D cases
    if (smeared_ && SmearObj !=NULL) attach_smearing(SmearObj);
  }

  ~Action_Nf(){}

  void init(const RandNum& rand);
  void observer_update();

  double calc_H();
  GaugeField md_force();

};
#endif
