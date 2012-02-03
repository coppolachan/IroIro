/*! @file  eoUtils.cpp
 *  @brief definition of utilities related to even/odd preconditioning
 */
#include "eoUtils.hpp"

namespace EvenOddUtils{

  void Inverter_WilsonLike::invert(Field& sol,const Field& src){

    Field be = D_->mult_ee_inv(Field(src[esub_]));
    Field bo = D_->mult_oo_inv(Field(src[osub_]));
    
    be -= D_->mult_eo(bo);

    Field ye(fsize_);
    SolverOutput monitor = slv_->solve(ye,D_->mult_dag(be));
#if VERBOSITY > 0
    monitor.print();
#endif
    bo -= D_->mult_oe(ye);
    
    for(int hs=0; hs<CommonPrms::instance()->Nvol()/2; ++hs){
      sol.set(ff_.islice(idx_eo_->esec(hs)), ye[ff_.islice(hs)]);
      sol.set(ff_.islice(idx_eo_->osec(hs)), bo[ff_.islice(hs)]);
    }
  } 
  
  void Inverter_WilsonLike::invert_dag(Field& sol,const Field& src){

  }
}
