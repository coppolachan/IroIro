/*! @file  eoUtils.cpp
 *  @brief definition of utilities related to even/odd preconditioning
 */
#include "eoUtils.hpp"
#include<cassert>

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

namespace EvenOddUtils{

  //Returns the multiplication of D*(full vector) in 5d
  Field Inverter_WilsonLike::mult(const Field& f) const{

    Field be = D_->mult_ee(Field(f[esub_]));
    be += D_->mult_ee(D_->mult_eo(Field(f[osub_])));//see def of mult_eo

    Field bo = D_->mult_oo(Field(f[osub_]));
    bo += D_->mult_oo(D_->mult_oe(Field(f[esub_]))); 

    Field out(ff_.size());
    // when SiteIndexEO is used this is equivalent to merging the even/odd slices of the 4d volumes in 
    // sequence like this (s index in the 5d slice if Nex=/=1)
    //
    // EvenVec[s=0 | s=1 | s=2 | ...] OddVec[s=0 | s=1 | s=2 | ...]
    // merged to
    // FullVec[ s=0 Even | s=0 Odd | s=1 Even | s=1 Odd | s=2 Even | s=2 Odd | ... ]
    // otherwise is using the full indexing (lexicographic) interleaving even and odd sites
#pragma omp parallel
    {
#pragma omp for
      for(int ex=0; ex<Nex_; ++ex){
	for(int hs=0; hs<Nvh_; ++hs){
	  out.set(ff_.islice(idx_->esec(hs),ex), be[fh_.islice(hs,ex)]);
	  out.set(ff_.islice(idx_->osec(hs),ex), bo[fh_.islice(hs,ex)]);
	}
      } 
    }

    return out;
  }

  void Inverter_WilsonLike::invert(Field& sol,const Field& src)const{
    // Inverting the 5d part of the 4d operator

    assert(sol.size()==ff_.size());
    assert(src.size()==ff_.size());
    
    long double timing;
    FINE_TIMING_START(timing);

    Field be = D_->mult_ee_inv(Field(src[esub_]));
    Field bo = D_->mult_oo_inv(Field(src[osub_]));

    be -= D_->mult_eo(bo);

    
    Field ye(ff_.size());
    SolverOutput monitor = slv_->solve(ye,D_->mult_dag(be));// to get y=M^-1 b
 
#if VERBOSITY > 0
    monitor.print();
#endif
    
    bo -= D_->mult_oe(ye);

    for(int ex=0; ex<Nex_; ++ex){
#pragma omp parallel
      {
#pragma omp for
	for(int hs=0; hs<Nvh_; ++hs){
	  sol.set(ff_.islice(idx_->esec(hs),ex), ye[fh_.islice(hs,ex)]);
	  sol.set(ff_.islice(idx_->osec(hs),ex), bo[fh_.islice(hs,ex)]);
	}
      }
      
    }

   FINE_TIMING_END(timing);
    _Message(TIMING_VERB_LEVEL,"[Timing] 4d DWF Inverter_WilsonLike solve :"<<timing<<"\n"); 
  }

  void Inverter_WilsonLike::test(Field& Df, const Field& f)const{
    
    Field Dfe = D_->mult_ee(Field(f[esub_]));
    Dfe += D_->mult_ee(D_->mult_eo(Field(f[osub_])));

    Field Dfo = D_->mult_oo(Field(f[osub_]));
    Dfo += D_->mult_oo(D_->mult_oe(Field(f[esub_])));

    for(int i=0;i<Dfe.size()/2;++i)
      CCIO::cout<<" Dfe["<<i<<"]=("<< Dfe[2*i]<<","<<Dfe[2*i+1]<<")"
		<<" Dfo["<<i<<"]=("<< Dfo[2*i]<<","<<Dfo[2*i+1]<<")"
		<<std::endl;

    for(int ex=0; ex<Nex_; ++ex){
      for(int hs=0; hs<Nvh_; ++hs){
	Df.set(ff_.islice(idx_->esec(hs),ex), Dfe[fh_.islice(hs,ex)]);
	Df.set(ff_.islice(idx_->osec(hs),ex), Dfo[fh_.islice(hs,ex)]);
      }
    } 
  }

}
