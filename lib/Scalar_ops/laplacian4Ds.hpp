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
#include "include/timings.hpp"
#include <omp.h>
#include <cassert>

template <typename FMT>
class Laplacian4Ds:public ScalarOp {
  const Field * const u_;
  std::vector<int> tsl_size_;
  Format::Format_G gf_;
  FMT ff_;
  void setup();

  inline int sr(int s,int c)const{return 2*(s*NC_+c);}
  inline int si(int s,int c)const{return 2*(s*NC_+c)+1;}
  inline int re(int c1,int c2)const{return 2*(c1*NC_+c2);}
  inline int im(int c1,int c2)const{return 2*(c1*NC_+c2)+1;}

public:
  Laplacian4Ds(Field* u)
  :u_(u),gf_(CommonPrms::instance()->Nvol()),
   ff_(CommonPrms::instance()->Nvol()){ setup();}

  const Field mult(const Field&)const;
  double func(double x)const{ return x;}
  size_t fsize()const {return ff_


.size();}
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
    GeneralField<Field,FMT,OneDimTag> tmp;
#ifdef IBM_BGQ_WILSON
    GeneralField<Field,FMT,OneDimTag> psi(F);
    GeneralField<Field,FMT,OneDimTag> tmp2(F);
#pragma omp parallel
    {
      shiftField(psi,F.data.getaddr(0),i, Forward());
    }
#else
    GeneralField<Field,FMT,OneDimTag> psi(shiftField(F,i,Forward()));
#endif
    

    
    for(int t=0; t<CommonPrms::Nt(); ++t){
#pragma omp parallel if(tsl_size_[t]%omp_get_num_threads() == 0)
      {
	int ns = tsl_size_[t]/omp_get_num_threads();
	int is = omp_get_thread_num()*ns;
	

	for(int n=is; n<is+ns; ++n){
	  int site = SiteIndex::instance()->slice_t(t,n);

	  //Gauge field (su3 18 doubles)
	  double *up = const_cast<Field*>(u_)->getaddr(gf_.index(0,site,i));
	  
	  // Fermion/Scalar field
	  double *fp = const_cast<Field&>(f).getaddr(ff_.index(0,site));

	  // Fermion/Scalar field
	  double *pp = psi.data.getaddr(ff_.index(0,site));
	  
	  double *c1p = chi1.data.getaddr(ff_.index(0,site));
	  double *ttp = tmp.data.getaddr(ff_.index(0,site));
	  
	  int Nd = ff_.Nin()/NC_/2;
	  
	  // Matrix vector product
#pragma disjoint(*up,*fp, *pp, *c1p, *ttp)

	  for(int s=0; s<Nd; ++s){
#pragma unrollandfuse(NC_)
	    for(int c=0; c<NC_; ++c){
	      int sc = sr(s,c);
	      for(int c1=0; c1<NC_; ++c1){
		int sri= sr(s,c1);
		
		c1p[sc]   += up[re(c,c1)]*pp[sri]   -up[im(c,c1)]*pp[sri+1];
		c1p[sc+1] += up[re(c,c1)]*pp[sri+1] +up[im(c,c1)]*pp[sri];
		
		ttp[sc]   += up[re(c1,c)]*fp[sri]   +up[im(c1,c)]*fp[sri+1];
		ttp[sc+1] += up[re(c1,c)]*fp[sri+1] -up[im(c1,c)]*fp[sri];
	      }
	 

	    }
	  }  
	}
      } // parallel section

    }// t loop 
  
#ifdef IBM_BGQ_WILSON
#pragma omp parallel
    {
      shiftField(tmp2,tmp.data.getaddr(0),i, Backward());
    }
    chi2 +=tmp2;
#else
    chi2 += shiftField(tmp,i,Backward());
#endif

  }

  F *= -6.0;
  F += chi1;
  F += chi2;
  
  return F.data;    
}

typedef Laplacian4Ds<Format::Format_F>  Laplacian4DF;
typedef Laplacian4Ds<Format::Format_S>  Laplacian4DS;

#endif
