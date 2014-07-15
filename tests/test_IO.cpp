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

#define NO_APPEND_MODE false
#define APPEND_MODE    true

#define DO_SU3_CHECK true
#define NO_SU3_CHECK false

int Test_IO::run(){
  CCIO::cout << "Test IO Starts ------------------------\n";
  CCIO::cout << "Running on "<<Communicator::instance()->size() << " nodes\n";
  
  ////////////////////////////////
  uint64_t offset = 0;
  std::string conf_type = "NERSC";
  ////////////////////////////////
  CCIO::cout << "Testing "<< conf_type << " writer/reader \n";

  Staples Staple;
  CCIO::cout << "Plaquette : " << Staple.plaquette(Gfield_) << std::endl;

  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.data, "temp_out", NO_APPEND_MODE, conf_type);

  CCIO::ReadFromDisk< Format::Format_G >(Gfield_.data, "temp_out", offset, conf_type, DO_SU3_CHECK );
  Staples Staple2;
  CCIO::cout << "Plaquette : " << Staple2.plaquette(Gfield_) << std::endl;

  return 0;
}
