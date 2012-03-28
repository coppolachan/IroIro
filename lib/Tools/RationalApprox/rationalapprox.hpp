/*!

 * @file rationalapprox.hpp

 * @brief Calculates the rational approximation for a given exponential function

 */

#ifndef RATIONAL_APPROX_HPP_
#define RATIONAL_APPROX_HPP_

#include <vector>
#include "include/pugi_interface.h"

/*! @brief Remez approximation parameters
  
  It is assumed that the function to approximate is of the form \$x^{a/b}\$
*/
struct RationalApprox_params{
  int numerator_deg;   /*!< @brief degree of the numerator polynomial */
  int denominator_deg; /*!< @brief degree of the numerator polynomial */

  int exponent_num;    /*!< @brief exponent numerator (\$a\$) */ 
  int exponent_den;    /*!< @brief exponent denominator (\$b\$) */ 

  int gmp_remez_precision;  /*! @brief precision used by GMP Remez algorithm */
  double lambda_low;           /*! @brief lower boundary of approximation interval */
  double lambda_high;          /*! @brief upper boundary of approximation interval */
};



/*! Calculates rational approximation valid in \$[ \epsilon_{min}, 1]\$ of the form: 
                                                            
  \f[ f(x)= r_0 + \sum_{i=0}^{i<N_{max}} \frac{r_i}{x+ p_i} \f]

  where \$r_i\$ and \$p_i\$ are the residuals and poles of the approximation.
*/                                                         
class RationalApprox {
private:
  RationalApprox_params Params;

  // Parameters for rational expansion (i.e. force, pseudofermions,...)
  int approximation_order;
  double min_epsilon;
  double RA_a0; /*!< @brief Rational Approximation Constant term */
  std::vector<double> RA_res; /*!< @brief Rational Approximation Residuals */
  std::vector<double> RA_pole;  /*!< @brief Rational Approximation Poles */

  RationalApprox(); // hide default constructor

  void fill();
public:
  /*! XML Constructor */
  RationalApprox(XML::node Approx_node);

  /*! Standard Constructor */
  RationalApprox(RationalApprox_params Par);

  std::vector<double> Residuals() { return RA_res; }
  std::vector<double> Poles() { return RA_pole; }
  double Const() { return RA_a0; }

};

#endif
