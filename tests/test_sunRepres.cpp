//------------------------------------------------------------------------
/*!
 * @file test_sunRepres.cpp
 *
 * @brief run() function for sunRep classes test
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <time.h>

#include "test_sunRepres.hpp"
#include "Tools/sunRepresentations.hpp"
#include "Tools/sunMatUtils.hpp"

int Test_sunRep::run(){
  using namespace SUNmatUtils;

  CCIO::cout << "Test sunRep Starts ------------------------\n";
  CCIO::cout << "Running on "<<Communicator::instance()->size() << " nodes\n";

  CCIO::cout << ":::::: Group SU("<<NC_<<")\n"; 
 
  SUNRep<3> Representations;

  for (int a = 0; a < 8; a++)
    Representations.lambda_fund(a).print();
  
  
  for (int a = 0; a < 8; a++)
    Representations.lambda_adj(a).print();
  

  return 0;
}
