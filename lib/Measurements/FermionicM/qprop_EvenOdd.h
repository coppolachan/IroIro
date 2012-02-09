/*!
 * @file qprop_EvenOdd.h
 * @brief Declaration of Quark propagator class Qprop_EvenOdd
 */
#ifndef QPROP_EVENODD_INCLUDED
#define QPROP_EVENODD_INCLUDED

#include "quark_propagators.hpp"
#include "Dirac_ops/eoUtils.hpp"

class DiracWilsonLike_EvenOdd;
class Solver;
class Source;
/*!
 * @class Qprop_EvenOdd
 * @brief Calculates the Quark propagator \f$D^{-1}\f$ with even/odd precond.
 */
class Qprop_EvenOdd : public QuarkPropagator{
private:
  const EvenOddUtils::Inverter_WilsonLike InvD_; /*!< @brief it contains e/o solver */
  const int fsize_;
  int Nc_,Nd_;
public:
  Qprop_EvenOdd(const DiracWilsonLike_EvenOdd* D, const Solver* solver)
    :InvD_(D,solver),fsize_(2*D->fsize()),
     Nc_(CommonPrms::instance()->Nc()),
     Nd_(CommonPrms::instance()->Nd()){
    CCIO::cout<<"Qprop_EvenOdd created"<<std::endl;
}

  ~Qprop_EvenOdd(){}
  
  /*! @brief Calculates the propagator 
   * N.B. A typedef std::vector<Field> prop_t;
   * is defined in quark_propagators.hpp
   */
  void calc(prop_t& xq,Source& src) const;
  void calc(prop_t& xq,Source& src, int, int) const;
  void calc(prop_t& xq,const prop_t& prp) const;
};
#endif
