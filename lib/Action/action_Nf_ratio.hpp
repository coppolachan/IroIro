/*!
 * @file action_Nf_ratio.hpp
 * @brief Declaration of Action_Nf_ratio class
 */
#ifndef ACTION_NF_RATIO_INCLUDED
#define ACTION_NF_RATIO_INCLUDED

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
struct Action_Nf_ratio_params {
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
  
  Action_Nf_ratio_params(){};

  Action_Nf_ratio_params(const XML::node node)
    :degree_(3),precision_(3),b_low_(3),b_high_(3){
    
    XML::read(node, "Flavors", n_flav_, MANDATORY); 
    XML::read(node, "Pseudofermions", n_pseudof_, MANDATORY);
    
    XML::read_array(node, "ApproxDegree", degree_, MANDATORY);
    XML::read_array(node, "Precision", precision_, MANDATORY);
    XML::read_array(node, "BoundaryLow", b_low_, MANDATORY);
    XML::read_array(node, "BoundaryHigh", b_high_, MANDATORY);
  }
};




/*!
 * @class Action_Nf_ratio
 * @brief Class to calculate HMC action term \f$det(D1)/det(D2)\f$
 */
class Action_Nf_ratio : public Action{
private:
  GaugeField* const u_; /*!< @brief The gauge field */
  DiracWilsonLike* D1_; /*!< @brief The numerator kernel */
  DiracWilsonLike* D2_; /*!< @brief The denominator kernel */
  RationalSolver* slv1_;  
  RationalSolver* slv2_;
  Action_Nf_ratio_params Params_;
  const size_t fermion_size_;
  std::vector<Field> phi_; /*!< @brief Vector of pseudofermion fields */

  bool smeared_;
  SmartConf* smart_conf_;

  // Rational approximations numerator
  RationalApprox MetropolisApprox_;  
  RationalApprox MolecularDynApprox_;
  RationalApprox PseudoFermionsApprox_;

  // Rational approximations denominator
  RationalApprox MetropolisApprox_Den_;  
  RationalApprox MolecularDynApprox_Den_;

  Field DdagD1_inv(const Field& src);
  Field DdagD2_inv(const Field& src);

  void attach_smearing(SmartConf*);
public:
  Action_Nf_ratio(GaugeField* const GField, 
		  DiracWilsonLike* const D1,
		  DiracWilsonLike* const D2,
		  RationalSolver* Solv1,
		  RationalSolver* Solv2,
		  Action_Nf_ratio_params Par,
		  bool smeared = false,
		  SmartConf* smart_conf = NULL)
    :u_(GField),
     D1_(D1), 
     D2_(D2),
     slv1_(Solv1), 
     slv2_(Solv2),
     Params_(Par),
     smeared_(smeared),
     fermion_size_(D1->fsize()),
     phi_(Params_.n_pseudof_),
     smart_conf_(smart_conf),
     MetropolisApprox_(RationalApprox_params(Params_.degree_[MetroStep],
					     Params_.degree_[MetroStep],
					     Params_.n_flav_,
					     2*Params_.n_pseudof_,
					     Params_.precision_[MetroStep],
					     Params_.b_low_[MetroStep],
					     Params_.b_high_[MetroStep])),
     MetropolisApprox_Den_(RationalApprox_params(Params_.degree_[MetroStep],
						 Params_.degree_[MetroStep],
						 Params_.n_flav_,
						 4*Params_.n_pseudof_,
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
     MolecularDynApprox_Den_(RationalApprox_params(Params_.degree_[MDStep],
						   Params_.degree_[MDStep],
						   Params_.n_flav_,
						   4*Params_.n_pseudof_,
						   Params_.precision_[MDStep],
						   Params_.b_low_[MDStep],
						   Params_.b_high_[MDStep])),
     PseudoFermionsApprox_(RationalApprox_params(Params_.degree_[PFStep],
						 Params_.degree_[PFStep],
						 Params_.n_flav_,
						 4*Params_.n_pseudof_,
						 Params_.precision_[PFStep],
						 Params_.b_low_[PFStep],
						 Params_.b_high_[PFStep]))
{
    for(int i=0; i<Params_.n_pseudof_; ++i)
      phi_[i].resize(fermion_size_); //takes care of EvenOdd and 5D cases


    //Temporary
    if(smart_conf_!= NULL) 
      assert(u_== smart_conf_->get_current_conf());

    //It is assumed that the rational approximation is always
    //applied to (M^dag M)
    //See the factor 2 in the denominators
  }
  
  ~Action_Nf_ratio(){}
  
  void init(const RandNum& rand);  
  void observer_update();

  double calc_H();
  GaugeField md_force();
};

#endif
