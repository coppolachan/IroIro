/*! @file  eoUtils.hpp
 *  @brief utilities related to even/odd preconditioning
 */
#ifndef EOUTILS_INCLUDED
#define EOUTILS_INCLUDED

#include "Dirac_ops/dirac.hpp"
#include "Solver/solver.hpp"
#include "Main/Geometry/siteIndex_EvenOdd.hpp"
#include "include/format_F.h"
#include <valarray>

namespace EvenOddUtils{

  class Inverter_WilsonLike{
  private:
    const DiracWilsonLike_EvenOdd* D_;
    const Solver* slv_;
    const SiteIndex_EvenOdd* idx_;
    const Format::Format_F fh_;  /*!< @brief half-size format */
    const Format::Format_F ff_;  /*!< @brief full-size format */
    const std::valarray<size_t> esub_;/*!< @brief generalized slice for even sites*/
    const std::valarray<size_t> osub_;/*!< @brief generalized slice for odd sites*/
  public:
    Inverter_WilsonLike(const DiracWilsonLike_EvenOdd* D,const Solver* Solver)
      :D_(D),slv_(Solver),
       fh_(D->get_fermionFormat()),
       ff_(2*fh_.Nvol(),fh_.Nex()),
       idx_(SiteIndex_EvenOdd::instance()),
       esub_(ff_.get_sub(idx_->esec())),
       osub_(ff_.get_sub(idx_->osec())){}

    void invert(Field& sol,const Field& src) const;
    void test(Field& Df, const Field& f) const;
  };
}

#endif
