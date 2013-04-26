/*
  @file stoutSmear.hpp
  @brief Declares Stout smearing class
*/
#ifndef STOUT_SMEAR_H_
#define STOUT_SMEAR_H_

#include "include/commonPrms.hpp"
#include "include/common_fields.hpp"
#include "APEsmear.hpp"
#include <valarray>

/*!  @brief Stout smearing of link variable. */
class Smear_Stout: public Smear {
private:
  const std::valarray<double> d_rho;
  const Smear* SmearBase;

  double func_xi0(double w) const;
public:
  Smear_Stout(Smear* base):SmearBase(base){}

  /*! Default constructor */
  Smear_Stout():SmearBase(new Smear_APE()){}

  ~Smear_Stout(){}

  void smear(GaugeField&,const GaugeField&) const;
  void BaseSmear(GaugeField&, const GaugeField&) const;
  void derivative(GaugeField&, const GaugeField&, const GaugeField&) const;
  void exponentiate_iQ(GaugeField&, const GaugeField&) const;
};

#endif  
