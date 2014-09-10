#include "include/messages_macros.hpp"
#include "test_LapSource.hpp"
#include "lib/Fopr/fopr_Exp.h"
#include "common_fields.hpp"
#include "lib/Measurements/FermionicM/source_types.hpp"
#include "lib/Scalar_ops/laplacian4Ds.hpp"
#include "format_F.h"
#include <iostream>

using namespace std;

int Test_LapSource::run(){

  InputConfig config = input_.getConfig();
  Laplacian4Ds<Format::Format_F> Lap(config.getGconf());
  Fopr_Scalar Flap(&Lap);

  vector<int> sp;
  for(int d=0; d<CommonPrms::Ndim(); ++d) sp.push_back(0);
  int Nvol = CommonPrms::Nvol();
  Source_local<Format::Format_F> src(sp,Nvol);

  Format::Format_F ff(Nvol);

  int Ly = CommonPrms::Ly();
  int Lz = CommonPrms::Lz();

  //  int Nd = CommonPrms::Nd();
  //  int Nc = CommonPrms::Nc();
  //  for(int s=0; s<Nd; ++s){
  //    for(int c=0; s<Nc; ++c){
  
  int s=0, c=0;
  Field v =src.mksrc(s,c);

  for(int N=0; N<210; N+=10) {
    double alpha = 0.1*double(N);
      
    Fopr_Lexp eSmr(N,alpha,&Flap);
    Field smrd = eSmr.mult(v);            
    Field Lv = Lap.mult(v);
    
    CCIO::cout<<"N="<<N<<" alpha="<<alpha<<" smrd.norm()="<<smrd.norm() <<"\n";

    if(SiteIndex::instance()->global_t(0) == 0){
      for(int sl= 0; sl<SiteIndex::instance()->slsize(0,TDIR); ++sl){
	
	int site = SiteIndex::instance()->slice_t(0,sl);
	int gsite = SiteIndex::instance()->get_gsite(site);
	
	int gx = SiteIndex::instance()->g_x(gsite);
	int gy = SiteIndex::instance()->g_y(gsite);
	int gz = SiteIndex::instance()->g_z(gsite);

	/*
	double re = Lv[ff.index_r(c,s,site)];
	double im = Lv[ff.index_i(c,s,site)];
	std::cout<<gx<<" "<<gy<<" "<<gz<<"  "<<re<<"\n";
	*/

	double re = smrd[ff.index_r(c,s,site)];
	double im = smrd[ff.index_i(c,s,site)];
	Communicator::instance()->sync();
	std::cout<<gx<<" "<<gy<<" "<<gz<<"  "<< sqrt(re*re+im*im)<<"\n";
	Communicator::instance()->sync();
      }
    }
  }
}

