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

 public:
  Smear_APE(const std::valarray<double>& rho_, 
	    const Format::Format_G& gf):
    Ndim(CommonPrms::Ndim()),
    Gformat(gf),
    rho(rho_){}


  ~Smear_APE(){};

  //void set_weights(const std::valarray<double>&);
  void smear(Field&, const Field&);

};

#endif  
