//------------------------------------------------------------------------
/*!
 * @file test_IO.hpp
 *
 * @brief run() function for  input/output functions test
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------

#include <time.h>

#include "test_IO.hpp"
#include "include/format_G.h"
#include "Communicator/comm_io.hpp"
#include "Communicator/fields_io.hpp"
#include "Measurements/GaugeM/staples.h"

int Test_IO::run(XML::node node){
  CCIO::cout << "Test IO Starts ------------------------\n";
  CCIO::cout << "Running on "<<Communicator::instance()->size() << " nodes\n";
  
  Staples Staple(Gfield_.Format);
  CCIO::cout << "Plaquette : " << Staple.plaquette(Gfield_.U) << std::endl;

  //CCIO::SaveOnDisk< Format::Format_G >(Gfield_.U, "temp_out");



}
