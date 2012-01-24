#ifndef APE_SMEAR_H_
#define APE_SMEAR_H_

#include <valarray>

#include "BaseSmear.hpp"
#include "include/commonPrms.h"
#include "include/format_G.h"

class Field;

/*!
  @brief APE type smearing of link variables.
 */
class Smear_APE: public Smear{

 private:
  const int Ndim;
  const Format::Format_G& Gformat;
  const std::valarray<double> rho;/*!< Array of weights */

  //This member must be private - we do not want to control from outside 
  std::valarray<double> set_rho(const double)const ;

 public:
  Smear_APE(const std::valarray<double>& rho_, 
	    const Format::Format_G& gf):
    Ndim(CommonPrms::Ndim()),
    Gformat(gf),
    rho(rho_){}

  Smear_APE(double rho_val, 
	    const Format::Format_G& gf):
    Ndim(CommonPrms::Ndim()),
    Gformat(gf),
    rho(set_rho(rho_val)){}

  Smear_APE(const Format::Format_G& gf):
    Ndim(CommonPrms::Ndim()),
    Gformat(gf),
    rho(set_rho(1.0)){}

   ~Smear_APE(){};

  void smear(Field&, const Field&) const;
  void derivative(Field&,const Field&, const Field&) const;
};

#endif  
