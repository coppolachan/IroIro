//------------------------------------------------------------------------  
// fopr_eo.h       
//------------------------------------------------------------------------
#ifndef FOPR_EO_INCLUDED
#define FOPR_EO_INCLUDED

#ifndef FOPR_INCLUDED
#include "fopr.h"
#endif

template<typename EO_DIRAC> class Fopr_Dhat :public Fopr{
private:
  const EO_DIRAC* D_;
public:
  explicit Fopr_Dhat(const EO_DIRAC* D):D_(D){}

  const Field mult(const Field& f) const{
    Field w = D_->mult_ee(f);
    w -= D_->mult_eo(D_->mult_ooinv(D_->mult_oe(f)));
    return w;
  }
  const Field mult_dag(const Field& f) const{
    return D_->gamma5(mult(D_->gamma5(f)));
  }
  size_t fsize()const {return D_->fsize();}
};

template<typename EO_DIRAC> class Fopr_Dhat_dag :public Fopr{
private:
  const EO_DIRAC* D_;
public:
  explicit Fopr_Dhat_dag(const EO_DIRAC* D):D_(D){}

  const Field mult(const Field& f) const{
    return D_->gamma5(mult_dag(D_->gamma5(f)));
  }
  const Field mult_dag(const Field& f) const{
    Field w = D_->mult_ee(f);
    w -= D_->mult_eo(D_->mult_ooinv(D_->mult_oe(f)));
    return w;
  }
  size_t fsize()const {return D_->fsize();}
};

template<typename EO_DIRAC> class Fopr_Hhat :public Fopr{
private:
  const EO_DIRAC* D_;
public:
  explicit Fopr_Hhat(const EO_DIRAC* D):D_(D){}

  const Field mult(const Field& f) const{
    Field w = D_->mult_ee(f);
    w -= D_->mult_eo(D_->mult_ooinv(D_->mult_oe(f)));
    return D_->gamma5(w);
  }
  const Field mult_dag(const Field& f) const{return mult(f);}
  size_t fsize()const {return D_->fsize();}
};
#endif
