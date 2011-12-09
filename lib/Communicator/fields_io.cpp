/*! 
 * @file fields_io.cpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#include "fields_io.hpp" 
#include "Main/Geometry/siteIndex.h"

namespace CCIO
{
  int SaveOnDisk(const std::vector<Field>& f, const char* filename) {
    std::ofstream Outputfile;

    if(Communicator::instance()->primaryNode()) {
      Outputfile.open(filename, std::fstream::out | std::fstream::binary);
      
      if (Outputfile.good()) {
	
	for (int i = 0; i < f.size(); ++i) {
	  for (int idx = 0; idx < f[i].getva().size(); ++idx){
	    Outputfile << f[i].getva()[idx];
	    
	  }
	}
	
	Outputfile.close();
      }

    }
    return 0;
  }

  template <typename T>
  int SaveOnDisk(const Field& f, const char* filename) {
    std::ofstream Outputfile;
    T format(CommonPrms::instance()->Nvol());
    Communicator* comm = Communicator::instance();
    size_t local_index, global_index;
    comm->sync();

    
    if(comm->primaryNode()) {
      Outputfile.open(filename, std::ios::binary);
      if (Outputfile.good()) {
	//Gather information
	Field Global(f.size()*CommonPrms::instance()->NP());
	std::valarray<double> local(f.size());
	//copy one by one the fields in the nodes into primary node
	for (int node=0; node < comm->size(); ++node) {
	  comm->send_1to1(local, f.getva(), f.size(), 0, node, node);
	  for (int site=0; site < format.Nvol(); ++site) {
	    for (int internal=0; internal< format.Nin(); ++internal) {
	       for (int external=0; external< format.Nex(); ++external) {
		 local_index = format.index(internal, site, external);
		 global_index = SiteIndex::instance()->
		   gsite(local_index);
		 
		 Global.set(global_index, f[local_index]);
	       }
	    }
	  }
	}// closing loop among processing nodes
	std::cout << "Binary writing "<<f.size()*sizeof(double)
		  <<" bytes on "<< filename << "\n";
	Global.write_stream(Outputfile);
      }
      Outputfile.close();
    }
    
    return 0;
  }
   
}
