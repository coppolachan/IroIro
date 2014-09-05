/*!
  @file laplacian4Ds.hpp
  @brief declares Laplacian4Ds class, which deals with 
  the application of Laplacian to each time-slice of 4D-spinors.
  As the fermion format FMT, Format_F or Format_S is assumed.
*/
#ifndef LAPLACIAN4DS_INCLUDED
#define LAPLACIAN4DS_INCLUDED

#include "scalarOp.hpp"
#include "commonPrms.hpp"
#include "format_G.h"
#include "Geometry/shiftField.hpp"
#include <cassert>

template <typename FMT>
class Laplacian4Ds:public ScalarOp {
  const Field * const u_;
  std::vector<int> tsl_size_;
  Format::Format_G gf_;
  FMT ff_;
  void setup();

  int sr(int s,int c)const{return 2*(s*NC_+c);}
  int si(int s,int c)const{return 2*(s*NC_+c)+1;}
  int re(int c1,int c2)const{return 2*(c1*NC_+c2);}
  int im(int c1,int c2)const{return 2*(c1*NC_+c2)+1;}

public:
  Laplacian4Ds(Field* u)
  :u_(u),gf_(CommonPrms::instance()->Nvol()),
   ff_(CommonPrms::instance()->Nvol()){ setup();}

  const Field mult(const Field&)const;
  double func(double x)const{ return x;}
  size_t fsize()const {return ff_.size();}
};

template <typename FMT>
void Laplacian4Ds<FMT>::setup(){
  for(int t=0; t<CommonPrms::Nt(); ++t)
    tsl_size_.push_back(SiteIndex::instance()->slsize(t,TDIR));
}

template <typename FMT>
const Field Laplacian4Ds<FMT>::mult(const Field& f)const{
  using namespace Mapping;

  GeneralField<Field,FMT,OneDimTag> F(f); 
  GeneralField<Field,FMT,OneDimTag> chi1,chi2;
 
  for(int i=0; i<NDIM_-1; ++i){ 
    GeneralField<Field,FMT,OneDimTag> psi(shiftField(F,i,Forward()));
    GeneralField<Field,FMT,OneDimTag> tmp;

    for(int t=0; t<CommonPrms::Nt(); ++t){
      for(int n=0; n<tsl_size_[t]; ++n){
	int site = SiteIndex::instance()->slice_t(t,n);

	double *fp = const_cast<Field&>(f).getaddr(ff_.index(0,site));
	double *up = const_cast<Field*>(u_)->getaddr(gf_.index(0,site,i));
	double *pp = psi.data.getaddr(ff_.index(0,site));

	double *c1p = chi1.data.getaddr(ff_.index(0,site));
	double *ttp = tmp.data.getaddr(ff_.index(0,site));

	int Nd = ff_.Nin()/NC_/2;

	for(int s=0; s<Nd; ++s){
	  for(int c=0; c<NC_; ++c){
	    for(int c1=0; c1<NC_; ++c1){
	      c1p[sr(s,c)] += up[re(c,c1)]*pp[sr(s,c1)] -up[im(c,c1)]*pp[si(s,c1)];
	      c1p[si(s,c)] += up[re(c,c1)]*pp[si(s,c1)] +up[im(c,c1)]*pp[sr(s,c1)];

	      ttp[sr(s,c)] += up[re(c1,c)]*fp[sr(s,c1)] +up[im(c1,c)]*fp[si(s,c1)];
	      ttp[si(s,c)] += up[re(c1,c)]*fp[si(s,c1)] -up[im(c1,c)]*fp[sr(s,c1)];
	    }
	  }
	}   
      }
    }
    chi2 += shiftField(tmp,i,Backward());
  }

  F *= -6.0;
  F += chi1;
  F += chi2;

  return F.data;    
}

typedef Laplacian4Ds<Format::Format_F>  Laplacian4DF;
typedef Laplacian4Ds<Format::Format_S>  Laplacian4DS;

#endif
