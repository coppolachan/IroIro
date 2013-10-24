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
  int Nt = CommonPrms::instance()->Nt();
  int Nvol = CommonPrms::instance()->Nvol();

  FermionField1sp input(Nvol/Nt);
  CCIO::cout << "Site  =" <<SiteIndex::instance()->site(2,2,2,2) << "\n";
  CCIO::cout << "input size  =" << input.size() << "\n";
  input.data.set(input.format.index(0,SiteIndex::instance()->site(2,2,2,0)), 1.0); 

  for(int i=0; i<input.size(); ++i)
    CCIO::cout <<"in["<<i<<"] = "<< input.data[i] << "\n";
   
  for(int t=0; t<Nt; ++t){
    Laplacian LapH(t,&(Gfield_.data));
    Field output = LapH.mult(input.data);
    for(int i=0; i<output.size(); ++i)
      CCIO::cout <<"t="<<t<<" "<<"out["<<i<<"] = "<< output[i] << "\n";
  }
  return 0;
}
