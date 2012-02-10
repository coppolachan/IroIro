//------------------------------------------------------------------------
/*!
 * @file test_IO.cpp
 *
 * @brief run() function for  input/output functions test
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <time.h>

#include "test_IO.hpp"
#include "Measurements/GaugeM/staples.hpp"

int Test_IO::run(){
  CCIO::cout << "Test IO Starts ------------------------\n";
  CCIO::cout << "Running on "<<Communicator::instance()->size() << " nodes\n";
  
  Staples Staple(Gfield_.Format);
  CCIO::cout << "Plaquette : " << Staple.plaquette(Gfield_.U) << std::endl;

  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.U, "temp_out");

  CCIO::ReadFromDisk< Format::Format_G >(Gfield_.U, "temp_out");
  Staples Staple2(Gfield_.Format);
  CCIO::cout << "Plaquette : " << Staple2.plaquette(Gfield_.U) << std::endl;

  return 0;
}
