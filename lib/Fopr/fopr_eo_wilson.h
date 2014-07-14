//---------------------------------------------------------------------
// fopr_eo_wilson.h
//---------------------------------------------------------------------
#ifndef FOPR_EO_WILSON_INCLUDED
#define FOPR_EO_WILSON_INCLUDED

#ifndef FOPR_EO_INCLUDED
#include "fopr_eo.h"
#endif

#ifndef DIRAC_WILSON_INCLUDED
#include "dirac_wilson.h"
#endif

// specialization of the class template Fopr_Dhat 

template<> class Fopr_Dhat<Dirac_Wilson::EvenOdd> :public Fopr{
private:
  const Dirac_Wilson::EvenOdd* D_;
public:
  explicit Fopr_Dhat(const Dirac_Wilson::EvenOdd* D):D_(D){}

  const Field mult(const Field& f) const{
    Field w(f);
    w -= D_->mult_eo(D_->mult_oe(f));
    return w;
  }
  const Field mult_dag(const Field& f) const{
    return D_->gamma5(mult(D_->gamma5(f)));
  }
  size_t fsize()const {return D_->fsize();}
};

template<> class Fopr_Dhat_dag<Dirac_Wilson::EvenOdd> :public Fopr{
private:
  const Dirac_Wilson::EvenOdd* D_;
public:
  explicit Fopr_Dhat_dag(const Dirac_Wilson::EvenOdd* D):D_(D){}

  const Field mult(const Field& f) const{
    return D_->gamma5(mult(D_->gamma5(f)));
  }
  const Field mult_dag(const Field& f) const{
    Field w(f);
    w -= D_->mult_eo(D_->mult_oe(f));
    return w;
  }
  size_t fsize()const {return D_->fsize();}
};

template<> class Fopr_Hhat<Dirac_Wilson::EvenOdd> :public Fopr{
private:
  const Dirac_Wilson::EvenOdd* D_;
public:
  explicit Fopr_Hhat(const Dirac_Wilson::EvenOdd* D):D_(D){}

  const Field mult(const Field& f) const{
    Field w(f);
    w -= D_->mult_eo(D_->mult_oe(f));
    return D_->gamma5(w);
  }
  const Field mult_dag(const Field& f) const{return mult(f);}
  size_t fsize()const {return D_->fsize();}
};

#endif
