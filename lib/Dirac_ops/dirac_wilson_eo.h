/*
 @file dirac_wilson_eo.h
 
 @brief Definition of Even Odd wilson operator
*/

#ifndef DIRAC_WILSON_EO_INCLUDED
#define DIRAC_WILSON_EO_INCLUDED

#include "dirac_wilson.h"


struct Even_tag{};
struct Odd_tag{};

namespace Dw{
  class Even;
  class Odd;
}

class Dirac_Wilson_EvenOdd : public Dirac_Wilson {
protected:

  const Dw::Even* even_;
  const Dw::Odd* odd_;

  Dirac_Wilson_EvenOdd(double kpp,const Field* u,Even_tag)
    :kpp_(kpp),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     ff_(new ffmt_t(Nvol_)), 
     gf_(new gfmt_t(Nvol_)),
     sf_up_(new shift_eup(ff_)),
     sf_dn_(new shift_edn(ff_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()),
     even_(0),odd_(0){}

  Dirac_Wilson_EvenOdd(double kpp,const Field* u,Odd_tag)
    :kpp_(kpp),u_(u),
     Nvol_(CommonPrms::instance()->Nvol()/2),
     Ndim_(CommonPrms::instance()->Ndim()),
     ff_(new ffmt_t(Nvol_)), 
     gf_(new gfmt_t(Nvol_)),
     sf_up_(new shift_oup(ff_)),
     sf_dn_(new shift_odn(ff_)),
     fsize_(ff_->size()),
     gsize_(gf_->size()),
     even_(0),odd_(0){}

  ~Dirac_Wilson_EvenOdd(){
    delete ff_;
    delete gf_;
    delete sf_up_;
    delete sf_dn_;
    if(even_) delete even_;
    if(odd_)  delete odd_;
  }

  const Dw::Even* get_even();
  const Dw::Odd* get_odd();
  

}




namespace Dw{
  class Even:private Dirac_Wilson{
  private:
    SiteIndex_eo* idx_;
  public:
    Even(double mq,const Field* u):Dirac_Wilson(mq,u,Even_tag()),
				   idx_(SiteIndex_eo::instance()){}
    SUNmat u(int site,int dir) const{
      return SUNmat((*u_)[gf_->cslice(0,idx_->esec(site),dir)]);
    }
    SUNmat u_dag(int site,int dir) const{
      return SUNmat((*u_)[gf_->cslice(0,idx_->esec(site),dir)]).dag();
    }
    const Field mult_oe(const Field&) const;
    const Field mult_ee(const Field& f) const{return f;}
    const Field mult_eeinv(const Field& f) const{return f;}
  };
  //
  class Odd:private Dirac_Wilson{
  private:
    SiteIndex_eo* idx_;
  public:
    Odd(double mq,const Field* u):Dirac_Wilson(mq,u,Odd_tag()),
				  idx_(SiteIndex_eo::instance()){}
    SUNmat u(int site,int dir) const{
      return SUNmat((*u_)[gf_->cslice(0,idx_->osec(site),dir)]);
    }
    SUNmat u_dag(int site,int dir) const{
      return SUNmat((*u_)[gf_->cslice(0,idx_->osec(site),dir)]).dag();
    }
    const Field mult_eo(const Field&) const;
    const Field mult_oo(const Field& f) const{return f;}
    const Field mult_ooinv(const Field& f) const{return f;}
  };
}



#endif
