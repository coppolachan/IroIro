/*!
 * @file meson_correlator.hpp
 *
 *
 * @brief Class for calculation of generic meson correlator
 *
 */

#ifndef MESON_CORRELATOR_HPP_
#define MESON_CORRELATOR_HPP_



typedef std::vector<Field> prop_t;
 
class MesonCorrelator {
  int Nc_;
  int Nd_;  
public:
  MesonCorrelator():Nc_(CommonPrms::instance()->Nc()),
		    Nd_(CommonPrms::instance()->Nd()){}
  const std::vector<double> calculate(const prop_t&, const prop_t&);

}


#endif //MESON_CORRELATOR_HPP_
