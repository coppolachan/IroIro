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
  BinaryFileReader::BinaryFileReader(const std::string& filename){
    //Initialize variables
    NumProcesses_ = Communicator::instance()->size();
    NodeID_       = Communicator::instance()->nodeid();
    isPrimary_    = Communicator::instance()->primaryNode();
    local_lattice_.resize(CommonPrms::instance()->Ndim());
    global_lattice_.resize(CommonPrms::instance()->Ndim());


  }




  int SaveOnDisk(std::vector<Field>& f, const char* filename) {
    std::fstream Outputfile;

    if(Communicator::instance()->primaryNode()) {
      Outputfile.open(filename, std::ios::out | std::ios::binary);
      
      if (Outputfile.is_open()) {
	Outputfile.write(reinterpret_cast<char* >(&f[0]), sizeof(double)*f[0].size()*f.size());
	Outputfile.close();
      }

    }
    return 0;
  }

}
