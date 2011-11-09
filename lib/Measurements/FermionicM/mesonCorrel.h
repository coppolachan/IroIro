//----------------------------------------------------------------------
// mesonCorrel.h
//----------------------------------------------------------------------
#ifndef MESONCORREL_INCLUDED
#define MESONCORREL_INCLUDED

#ifndef COMMONPRMS_INCLUDED
#include "commonPrms.h"
#endif

#ifndef COMMUNICATOR_INCLUDED
#include "communicator.h"
#endif

#include <vector>

using CommunicatorItems::pprintf;

class SiteIndex;
class Field;
typedef std::vector<Field> prop_t;

template<typename FMT> class MesonCorrel{
private:
  const SiteIndex* idx_;
  const FMT* fmt_;
  int Nx_,Ny_,Nz_,Nt_;
  int Nc_;
  int Nd_;  
  int Nvol_;
public:
  MesonCorrel():idx_(SiteIndex::instance()),
		Nx_(CommonPrms::instance()->Nx()),
		Ny_(CommonPrms::instance()->Ny()),
		Nz_(CommonPrms::instance()->Nz()),
		Nt_(CommonPrms::instance()->Nt()),
		Nc_(CommonPrms::instance()->Nc()),
		Nd_(CommonPrms::instance()->Nd()),
		Nvol_(CommonPrms::instance()->Nvol()),
		fmt_(new FMT(Nvol_)){}

  ~MesonCorrel(){delete fmt_;}

  const std::vector<double>   pp(const prop_t&, const prop_t&);
  /*
  const std::vector<double>  a4p(const prop_t&, const prop_t&);
  const std::vector<double>  pa4(const prop_t&, const prop_t&);
  const std::vector<double> a4a4(const prop_t&, const prop_t&);
  const std::vector<double> v1v1(const prop_t&, const prop_t&);
  const std::vector<double> v2v2(const prop_t&, const prop_t&);
  const std::vector<double> v3v3(const prop_t&, const prop_t&);
  */
};  

template<typename FMT> 
const std::vector<double> 
MesonCorrel<FMT>::pp(const prop_t& q1, const prop_t& q2){
  pprintf("contraction to make up  pp-correlator\n");

  std::vector<double> correl_local(Nt_,0.0);
  for(int site=0;site<Nvol_;++site){
    int t = idx_->t(site);
    for(int s1=0;s1<Nd_;++s1){
      for(int c1=0;c1<Nc_;++c1){
	//    {int s1=0;
	//      { int c1=0;
	for(int s2=0;s2<Nd_;++s2){
	  correl_local[t]+=(q1[c1+Nc_*s1][fmt_->cslice(s2,site)]
			    *q2[c1+Nc_*s1][fmt_->cslice(s2,site)]).sum();
	}
      }
    }
  }

  int Lt = CommonPrms::instance()->Lt();
  int ipet = Communicator::instance()->ipe(3);

  std::vector<double> correl_tmp(Lt,0.0);
  for(int t = 0; t< Nt_; ++t) correl_tmp[t +ipet*Nt_] = correl_local[t];

  std::vector<double> correl(Lt);
  for(int t = 0; t< Lt; ++t) 
    correl[t] = Communicator::instance()->reduce_sum(correl_tmp[t]);
  
  return correl;
}

#endif

