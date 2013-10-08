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
#include "Laplacian/laplacian.hpp"


int Test_LapH::run(){
  CCIO::cout << "Starting LapH" << std::endl;
   
  Laplacian LapH(&Gfield_);
  
  int Nvol = CommonPrms::instance()->Nvol();
  
  FermionField1sp input;
  CCIO::cout << "Site  =" <<SiteIndex::instance()->site(2,2,2,2) << "\n";
  CCIO::cout << "input size  =" << input.size() << "\n";
  input.data.set(input.format.index_r(0,SiteIndex::instance()->site(2,2,2,2)), 1.0); 
  /*
  for (int i = 0; i < input.size(); i++){
    CCIO::cout <<"in["<<i<<"] = "<< input.data[i] << "\n";
  }
  */
 
  
  FermionField1sp output = LapH.apply(input);
  

  for (int i = 0; i < output.size(); i++){
    CCIO::cout <<"out["<<i<<"] = "<< output.data[i] << "\n";
  }


  return 0;
}
