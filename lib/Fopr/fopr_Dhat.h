//------------------------------------------------------------------------
// fopr_Hhat.h                                                 
//------------------------------------------------------------------------
#ifndef FOPR_DHAT_INCLUDED
#define FOPR_DHAT_INCLUDED

#ifndef FOPR_INCLUDED
#include "fopr.h"
#endif

template<typename DIRAC> class Fopr_Dhat:public Fopr{
private:
  const DIRAC* D_;
  const bool myD_;
public:
  Fopr_Hhat(const typename DIRAC::Prms& prms,const Field* const u)
    :D_(new DIRAC(prms,u)),myD_(true){}
  explicit Fopr_Dhat(const DIRAC* D):D_(D),myD_(false){}
  ~Fopr_Dhat(){if(myD_) delete D_;}

  const Field mult(    const Field& f)const{
    Field w = D_->mult_ee(f);
    w-= D_->mult_eo(D_->mult_ooi(D_->mult_oe(f)));
    return w;
  }
  const Field mult_dag(const Field& f)const{return mult(f);}
  size_t fsize()const {return D_->fsize()/2;}
};

#endif
