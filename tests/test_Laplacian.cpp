//------------------------------------------------------------------------
/*!
 * @file test_Laplacian.cpp
 *
 * @brief run() function for Laplacian class test
 *
 * @author Jun-Ichi Noaki
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>

#include "test_Laplacian.hpp"
#include "Scalar_ops/laplacian.hpp"

int Test_LapH::run(){
  CCIO::cout << "Starting LapH" << std::endl;

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();

  int Nvol = CommonPrms::instance()->Nvol();

  SiteIndex3d idx;

  FermionField1sp input(Nvol/Nt);
  
  int bp[3] = {0,0,0};

  CCIO::cout << "input size =" << input.size() <<"\n";  
  CCIO::cout << "Site = ("<<bp[XDIR]<<","<<bp[YDIR]<<","<<bp[ZDIR]<<")\n";

  int nid = Communicator::instance()->nodeid(bp[XDIR]/Nx,bp[YDIR]/Ny,
					     bp[ZDIR]/Nz,0);

  Communicator::instance()->sync();
  int loc=0;
  if(Communicator::instance()->nodeid() == nid){
    loc = idx.site(bp[XDIR]%Nx,bp[YDIR]%Ny,bp[ZDIR]%Nz);
    input.data.set(input.format.index(0,loc),1.0); 
  }

  for(int site=0; site<input.size()/input.format.Nin(); ++site){
    for(int in=0; in<input.format.Nin()/2; ++in){
      std::cout <<"input["
		<<idx.global_x(idx.c_x(site))<<","
		<<idx.global_y(idx.c_y(site))<<","
		<<idx.global_z(idx.c_z(site))<<"]:"
		<<" c= "<<in<<" ("
		<< input[input.format.index_r(in,site)] << ","
		<< input[input.format.index_i(in,site)] << ")\n";
    }
  }
  Communicator::instance()->sync();

  for(int t=0; t<Nt; ++t){
    Laplacian LapH(t,&(Gfield_.data));
    Field output = LapH.mult(input.data);
    CCIO::cout<<"t="<<t<<"\n";
    for(int site=0; site<input.size()/input.format.Nin(); ++site){
      for(int in=0; in<input.format.Nin()/2; ++in){
	std::cout <<"output["
		  <<idx.global_x(idx.c_x(site))<<","
		  <<idx.global_y(idx.c_y(site))<<","
		  <<idx.global_z(idx.c_z(site))<<"]:"
		  <<" c= "<<in<<" ("
		  << output[input.format.index_r(in,site)] << ","
		  << output[input.format.index_i(in,site)] << ")\n";
      }
    }
    Communicator::instance()->sync();
  }
  return 0;
}
