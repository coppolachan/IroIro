/*!
  @file APEsmear.hpp
  @brief Declaration of Smear_APE class for APE smearing
*/

#ifndef APE_SMEAR_H_
#define APE_SMEAR_H_

#include <valarray>
#include "baseSmear.hpp"
#include "include/commonPrms.hpp"

/*!  @brief APE type smearing of link variables. */
class Smear_APE: public Smear{
private:
  const std::vector<double> rho;/*!< Array of weights */

  //This member must be private - we do not want to control from outside 
  std::vector<double> set_rho(const double)const ;

public:
  Smear_APE(const std::vector<double>& rho_):rho(rho_){}
  Smear_APE(double rho_val):rho(set_rho(rho_val)){}
  Smear_APE():rho(set_rho(1.0)){}
  ~Smear_APE(){}

  void smear(GaugeField&, const GaugeField&)const;
  void derivative(GaugeField&,const GaugeField&,const GaugeField&)const;
};

#endif  
