/*! @file  eoUtils.hpp
 *  @brief utilities related to even/odd preconditioning
 */
#include "Dirac_ops/dirac.h"
#include "include/format_F.h"

namespace EvenOddUtils{

  class Inverter_WilsonLike{

    const DiracWilsonLike_EvenOdd& D_;
    const int fsize_;
    const Solver& slv_;
    const Format::Format_F fh_;  /*!< @brief half-size format */
    const Format::Format_F ff_;  /*!< @brief full-size format */
    std::valarray<size_t> esub_; /*!< @brief generalized slice for even indices*/
    std::valarray<size_t> osub_; /*!< @brief generalized slice for for indices*/
  private:
    Inverter_WilsonLike(const DiracWilsonLike_EvenOdd& D,const Solver& Solver)
      :D_(D),fsize_(2*D.fsize()),
       slv_(Solver),
       fh_(D.get_fermionFormat()),
       ff_(2*(fh_.Nvol()),fh_.Nex()),
       idx_eo_(SiteIndex_eo::instance()),
       esub_(hf_.get_sub(idx_eo_->esec())),
       osub_(hf_.get_sub(idx_eo_->osec())){}

    void invert(Field& sol,const Field& src);
    void invert_dag(Field& sol,const Field& src);
  };

}
