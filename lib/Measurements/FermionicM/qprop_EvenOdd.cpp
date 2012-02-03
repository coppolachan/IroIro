//------------------------------------------------------------------------
/*!
 * @file qprop_EvenOdd.cpp
 * @brief Definition of Qprop_EvenOdd member functions 
 */
//-------------------------------------------------------------------------
#include "qprop_EvenOdd.h"
#include "Measurements/FermionicM/source.hpp"
#include "include/format_F.h"
#include "Communicator/comm_io.hpp"

using namespace std;
using namespace Format;

void Qprop_EvenOdd::calc(prop_t& xq,Source& src) const{

  xq.clear();
  SolverOutput monitor;
  Field ye(fsize_);
  
  Format::Format_F ff(2*fsize_);  

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      CCIO::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;
      // following block can be global

      Field be = D_->mult_ee_inv(src.mksrc(idx_eo_->esec(),s,c));
      Field bo = D_->mult_oo_inv(src.mksrc(idx_eo_->osec(),s,c));
	
      be -= D_->mult_eo(bo);
	
      monitor = slv_->solve(ye,D_->mult_dag(be));
#if VERBOSITY > 0
      CCIO::cout<<"s="<<s<<" c="<<c<< std::endl;
      monitor.print();
#endif
      bo -= D_->mult_oe(ye);
      Field sol(2*fsize_);
	
      for(int hs=0; hs<CommonPrms::instance()->Nvol()/2; ++hs){
	sol.set(ff.islice(idx_eo_->esec(hs)), ye[ff.islice(hs)]);
	sol.set(ff.islice(idx_eo_->osec(hs)), bo[ff.islice(hs)]);
      }
      xq.push_back(sol);
    }
  }
}

void Qprop_EvenOdd::calc(prop_t& xq,Source& src, int Nd, int Nc) const{}

void Qprop_EvenOdd::calc(prop_t& xq,const prop_t& prp)const{

  xq.clear();
  SolverOutput monitor;
  Field ye(fsize_);
  
  Format::Format_F ff(CommonPrms::instance()->Nvol());
  Format::Format_F hf(CommonPrms::instance()->Nvol()/2);
  /*
  vector<int> esub = hf.get_sub(idx_eo_->esec());
  vector<int> osub = hf.get_sub(idx_eo_->osec());
  */
  valarray<size_t> esub = hf.get_sub(idx_eo_->esec());
  valarray<size_t> osub = hf.get_sub(idx_eo_->osec());

  for(int s=0; s<Nd_;++s){
    for(int c=0; c<Nc_;++c){
      CCIO::cout<<" Dirac index ="<<s<<" Color ="<<c<<std::endl;
       
      Field be = D_->mult_ee_inv(Field(prp[Nc_*s+c][esub]));
      Field bo = D_->mult_oo_inv(Field(prp[Nc_*s+c][osub]));

      be -= D_->mult_eo(bo);

      monitor = slv_->solve(ye,D_->mult_dag(be));
#if VERBOSITY > 0
      CCIO::cout<<"s="<<s<<" c="<<c<< std::endl;
      monitor.print();
#endif
      bo -= D_->mult_oe(ye);
      Field sol(2*fsize_);

      for(int hs=0; hs<fsize_;++hs){
	sol.set(ff.islice(idx_eo_->esec(hs)), ye[ff.islice(hs)]);
	sol.set(ff.islice(idx_eo_->osec(hs)), bo[ff.islice(hs)]);
      }
      xq.push_back(sol);
    }
  }
}
