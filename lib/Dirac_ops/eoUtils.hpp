/*! @file  eoUtils.hpp
 *  @brief utilities related to even/odd preconditioning
 Time-stamp: <2013-06-20 13:07:33 noaki>
 */
#ifndef EOUTILS_INCLUDED
#define EOUTILS_INCLUDED

#include "Dirac_ops/dirac_WilsonLike.hpp"
#include "Solver/solver.hpp"
#include "Geometry/siteIndex_EvenOdd.hpp"
#include "include/format_F.h"
#include <valarray>

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

namespace EvenOddUtils{

  class Inverter_WilsonLike{
  private:
    const DiracWilsonLike_EvenOdd* D_;
    const Solver* slv_;
    int Nvol_,Nvh_,Nex_;
    const ffmt_t ff_;  /*!< @brief full-size format */
    const ffmt_t fh_;  /*!< @brief half-size format */
    const SiteIndex_EvenOdd* idx_;
    const std::valarray<size_t> esub_;/*!< @brief generalized slice for even sites*/
    const std::valarray<size_t> osub_;/*!< @brief generalized slice for odd sites*/
  public:
    Inverter_WilsonLike(const DiracWilsonLike_EvenOdd* D,const Solver* Solver)
      :D_(D),slv_(Solver),
       Nvol_(CommonPrms::instance()->Nvol()),Nvh_(Nvol_/2),
       Nex_(D->fsize()/Nvh_/ffmt_t::Nin()),
       ff_(Nvol_,Nex_),
       fh_(Nvh_,Nex_),
       idx_(SiteIndex_EvenOdd::instance()),
       esub_(ff_.get_sub(idx_->esec())),
       osub_(ff_.get_sub(idx_->osec())){ assert(D_);}

    void invert(Field& sol,const Field& src) const;
    
    const Field* getGaugeField_ptr()const{ return D_->getGaugeField_ptr(); }
    Field mult(const Field& in) const;
    void test(Field& Df, const Field& f) const;
  };
}

#endif
