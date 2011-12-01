//------------------------------------------------------------------------
/*!
 *
 * @file qprop_EvenOdd.cpp
 *
 * @brief Definition of Qprop_EvenOdd member functions 
 *
 */
//-------------------------------------------------------------------------
#include "qprop_EvenOdd.h"
#include "Measurements/FermionicM/source.h"
#include "include/format_F.h"
#include "Communicator/comm_io.hpp"

using namespace std;
using namespace Format;

void Qprop_EvenOdd::calc(prop_t& xq,Source& src) const{

  xq.clear();
  int Nconv;
  double diff;
  Field ye(fsize_);
  
  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
 
      if(Communicator::instance()->nodeid()==0) 
	CCIO::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;

      Field be(src.mksrc(idx_eo_->esec(),s,c));
      Field bo(src.mksrc(idx_eo_->osec(),s,c));
      
      be -= D_->mult_eo(D_->mult_oo_inv(bo));

      slv_->solve(ye,D_->mult_dag(be),diff,Nconv);

      CCIO::cout<<"s="<<s<<" c="<<c
		<<" Nconv= "<<Nconv 
		<<" diff= "<<diff<<std::endl;     

      bo -= D_->mult_oo_inv(D_->mult_oe(ye));

      Field sol(2*fsize_);
      Format_F ff(2*fsize_);

      for(int hs=0; hs<CommonPrms::instance()->Nvol()/2; ++hs){
	sol.set(ff.islice(idx_eo_->esec(hs)), ye[ff.islice(hs)]);
	sol.set(ff.islice(idx_eo_->osec(hs)), bo[ff.islice(hs)]);
      }
      xq.push_back(sol);
    }
  }
}

void Qprop_EvenOdd::calc(prop_t& xq,Source& src, int Nd, int Nc) const{};

void Qprop_EvenOdd::calc(prop_t& xq,const prop_t& prp)const{

  xq.clear();
  int Nconv;
  double diff;
  Field ye(fsize_);

  Format_F ff(2*fsize_);
  Format_F hf(fsize_);
  /*
  vector<int> esub = hf.get_sub(idx_eo_->esec());
  vector<int> osub = hf.get_sub(idx_eo_->osec());
  */
  valarray<size_t> esub = hf.get_sub(idx_eo_->esec());
  valarray<size_t> osub = hf.get_sub(idx_eo_->osec());

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){

      /*      
      Field be(fsize_);
      Field bo(fsize_);
      for(int i=0;i<fsize_;++i){
	be.set(i, prp[Nc_*s+c][esub[i]]);
	bo.set(i, prp[Nc_*s+c][osub[i]]);
      }
      */
      Field be(prp[Nc_*s+c][esub]);
      Field bo(prp[Nc_*s+c][osub]);

      be -= D_->mult_eo(D_->mult_oo_inv(bo));

      slv_->solve(ye,D_->mult_dag(be),diff,Nconv);

      CCIO::cout<<"s="<<s<<" c="<<c
		<<" Nconv= "<<Nconv 
		<<" diff= "<<diff<<std::endl;     

      bo -= D_->mult_oo_inv(D_->mult_oe(ye));

      Field sol(2*fsize_);

      for(int hs=0; hs<fsize_;++hs){
	//CCIO::cout<<"hs="<<hs<<endl;
	sol.set(ff.islice(idx_eo_->esec(hs)), ye[ff.islice(hs)]);
	sol.set(ff.islice(idx_eo_->osec(hs)), bo[ff.islice(hs)]);
      }
      xq.push_back(sol);
    }
  }
}
