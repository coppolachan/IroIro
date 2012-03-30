/*! 
 * @file action_Nf2.h
 * @brief Declaration of Action_Nf class
 *
 * Any number of flavours
 */
#ifndef ACTION_NF_INCLUDED
#define ACTION_NF_INCLUDED

#include <vector>

#include "Action/action.hpp"
#include "Dirac_ops/dirac.hpp"
#include "Solver/rationalSolver.hpp"
#include "Smearing/smartConf.hpp"
#include "Tools/RationalApprox/rationalapprox.hpp"


enum {MetroStep, MDStep, PFStep};

/*!
 * @brief Parameter container for Action_Nf class
 */
struct Action_Nf_params {
  int n_pseudof_;        /*!< @brief Number of pseudofermion fields */
  int n_flav_;           /*!< @brief Number of flavors */ 
  std::vector<int> degree_;   /*!< @brief Polynomial degree of approximation -
			   3 integers for Metropolis, Molecular Dynamics
			   and Pseudofermion steps, respectively */
  std::vector<int> precision_;/*!< @brief Precision for GMP routines -
			   3 integers for Metropolis, Molecular Dynamics
			   and Pseudofermion steps, respectively */
  std::vector<double> b_low_; /*!< @brief Lower boundary of approximation
			   interval -
			   3 integers for Metropolis, Molecular Dynamics
			   and Pseudofermion steps, respectively */
  std::vector<double> b_high_;/*!< @brief UpperLower boundary of approximation
			   interval -
			   3 integers for Metropolis, Molecular Dynamics
			   and Pseudofermion steps, respectively */
  
  Action_Nf_params(const XML::node node)
    :degree_(3),precision_(3),b_low_(3),b_high_(3){
    
    XML::read(node, "Flavors", n_flav_, MANDATORY); 
    XML::read(node, "Pseudofermions", n_pseudof_, MANDATORY);
    
    XML::read_array(node, "ApproxDegree", precision_, MANDATORY);
    XML::read_array(node, "BoundaryLow", b_low_, MANDATORY);
    XML::read_array(node, "BoundaryHigh", b_high_, MANDATORY);
  }
};


/*!
 * @class Action_Nf
 * @brief Class to calculate N-Flavors RHMC action term
 */
class Action_Nf :public Action{
private:
  GaugeField* const u_;      /*!< The gauge field */
  DiracWilsonLike* const D_; /*!< Dirac Kernel operator */ 
  RationalSolver* slv_;        /*!< RationalSolver operator */
  std::vector<Field> phi_;
  Action_Nf_params Params_;
  int fermion_size_;
  bool smeared_;
  SmartConf* SmartField_;
  
  // Rational approximations
  RationalApprox MetropolisApprox_;  
  RationalApprox MolecularDynApprox_;
  RationalApprox PseudoFermionsApprox_;
  
  Field DdagD_inv(const Field& src);
  void attach_smearing(SmartConf*);
public:
  /*!
   * @brief Standard constructor 
   * CG solver is assumed
   */
  Action_Nf(GaugeField* const GField,
	    DiracWilsonLike* const D, 
	    RationalSolver* Solv,
	    Action_Nf_params Par,
	    bool smeared = false,
	    SmartConf* SmearObj = NULL)
    :u_(GField),
     D_(D),
     slv_(Solv),
     Params_(Par),
     smeared_(smeared),
     fermion_size_(D->fsize()),
     phi_(Params_.n_pseudof_),
     MetropolisApprox_(RationalApprox_params(Params_.degree_[MetroStep],
					     Params_.degree_[MetroStep],
					     Params_.n_flav_,
					     2*Params_.n_pseudof_,
					     Params_.precision_[MetroStep],
					     Params_.b_low_[MetroStep],
					     Params_.b_high_[MetroStep])),
     MolecularDynApprox_(RationalApprox_params(Params_.degree_[MDStep],
					       Params_.degree_[MDStep],
					       Params_.n_flav_,
					       2*Params_.n_pseudof_,
					       Params_.precision_[MDStep],
					       Params_.b_low_[MDStep],
					       Params_.b_high_[MDStep])),
     PseudoFermionsApprox_(RationalApprox_params(Params_.degree_[PFStep],
						 Params_.degree_[PFStep],
						 Params_.n_flav_,
						 4*Params_.n_pseudof_,
						 Params_.precision_[PFStep],
						 Params_.b_low_[PFStep],
						 Params_.b_high_[PFStep])){
    for(int i=0; i<Params_.n_pseudof_; ++i)
      phi_[i].resize(fermion_size_); //takes care of EvenOdd and 5D cases
    if (smeared_ && SmearObj !=NULL) attach_smearing(SmearObj);
    
    //It is assumed that the rational approximation is always
    //applied to (M^dag M)
    //See the factor 2 in the denominators
  }

  ~Action_Nf(){}

  void init(const RandNum& rand);
  void observer_update();

  double calc_H();
  GaugeField md_force();

};
#endif
