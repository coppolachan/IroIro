/*! 
 * @file fields_io.cpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#include "fields_io.hpp" 

namespace CCIO
{
  int SaveOnDisk(const std::vector<Field>& f, const char* filename) {
    std::fstream Outputfile;

    if(Communicator::instance()->primaryNode()) {
      Outputfile.open(filename, std::fstream::out | std::fstream::binary);
      
      if (Outputfile.good()) {
	
	for (int i = 0; i < f.size(); i++) {
	  for (int idx = 0; idx < f[i].getva().size(); idx++){
	    Outputfile << f[i].getva()[idx];
	    
	  }
	}
	
	Outputfile.close();
      }

    }
    return 0;
  }
}
