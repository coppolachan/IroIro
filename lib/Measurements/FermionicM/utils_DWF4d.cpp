/*!
 * @file utils_DWF4d.cpp
 * @brief Declaration of utility functions for 4D-DWF op.
*/
#include "Measurements/FermionicM/utils_DWF4d.hpp"
#include "Dirac_ops/dirac_WilsonLike.hpp"
#include "Fields/field_expressions.hpp"

namespace Utils_DWF4D{

  /*!
   * @brief Calculates the sign function of the Kernel
   * It uses the equation
   * \f[{\rm sign}(H_{kernel}) = 2 \gamma_5 D(m=0) - \gamma_5\f]
   */
  const Field signKernel(const Dirac_DomainWall_4D& Ddwf,
			 const Field& f){
    using namespace FieldExpression;

    double mass = Ddwf.getMass();
    Field signK = Ddwf.mult(f);
    signK -= mass*f;
    signK /= 1.0 -mass; // (D(m) - m)/(1-m) = D(0)

    signK -= 0.5*f;
    signK = 2.0*Ddwf.gamma5(signK); // 2*\gamma_5(D(0)-1/2)
    return signK;
  }

  /*! @brief  Delta function = 1/4 * (1- sign^2(Hw)) */
  const Field delta(const Dirac_DomainWall_4D& Ddwf,const Field& phi){

    Field delta = signKernel(Ddwf,signKernel(Ddwf,phi));  //sign^2(Hw)
    delta -= phi;  //sign^2(Hw) -1  
    delta *= -0.25; // 1/4*(1-sign^2(Hw))
    return delta;
  }

}
