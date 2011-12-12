/*! 
 * @file fields_io.hpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#include "include/field.h"
#include "Main/Geometry/siteIndex.h"

#define CCIO_FILE_APPEND_MODE true

namespace CCIO {
  template <typename T>
  int SaveOnDisk(const Field& f, const char* filename, bool append_mode = false) 
  {
    T format(CommonPrms::instance()->Nvol());
    Field Global;
    Communicator* comm = Communicator::instance();
    size_t local_index, global_index;
    int num_nodes, total_volume;
    int* global_site = new int[format.Nvol()];
    
    total_volume = CommonPrms::instance()->Lvol();
    num_nodes    = CommonPrms::instance()->NP();
    
    //Initializations
    std::valarray<double> local(f.size());
    
    if(comm->primaryNode()) 
      Global.resize(f.size()*num_nodes);
    
    
    //Gather information
    for (int node=0; node < num_nodes; ++node) {
	comm->sync();
	//copy one by one the local fields stored in the nodes into primary node
	comm->send_1to1(local, f.getva(), local.size(), 0, node, node);
	
	for (int site = 0; site < format.Nvol(); ++site) {
	  global_site[site] = SiteIndex::instance()->
	    gsite(site);
	}
	comm->sync();    
	//copy global index array on primary node
	comm->broadcast(format.Nvol(), global_site[0], node);
	if(comm->primaryNode())
	  {	  
	    for (int site=0; site < format.Nvol(); ++site) {
	      for (int internal =0; internal < format.Nin(); ++internal) {
		for (int external =0; external < format.Nex(); ++external) {
		  local_index = format.index(internal, site, external);
		  //lime format: nb  external and site index are reverted from original order
		  global_index = internal +format.Nin()*(external+ format.Nex()*global_site[site]);
		  
		  Global.set(global_index, local[local_index]);
		}
	      }
	    }
	  }
      }
	   // closing loop among processing nodes
	   
   if(comm->primaryNode()) {
     std::ofstream Outputfile;
     if(append_mode) 
       Outputfile.open(filename, std::ios::binary | std::ios::app);
     else
       Outputfile.open(filename, std::ios::binary);
     
     if (Outputfile.good()) {
       
       std::cout << "Binary writing "<<f.size()*comm->size()*sizeof(double)
		 <<" bytes on "<< filename << "\n";
       Global.write_stream(Outputfile);
     }
     Outputfile.close();
   };
 
 delete[] global_site;
 
 return 0;
}
    
    template <typename T>
      int SaveOnDisk(const std::vector<Field>& f, const char* filename) {
      for (int field_num = 0; field_num < f.size(); ++field_num) {
	SaveOnDisk<T>(f[field_num], filename, CCIO_FILE_APPEND_MODE);
      }
    }


  template <typename T>
    int ReadFromDisk(Field& f, const char* filename, int offset = 0) {
    T format(CommonPrms::instance()->Nvol());
    Field Global;
    Communicator* comm = Communicator::instance();
    size_t local_index, global_index;
    int num_nodes, total_volume, local_volume;

    
    total_volume = CommonPrms::instance()->Lvol();
    local_volume = CommonPrms::instance()->Nvol();
    num_nodes    = CommonPrms::instance()->NP();
    
    //Read data from single file on disk
    if(comm->primaryNode()) {
      std::ifstream Inputfile;
      Global.resize(format.size()*num_nodes);

      Inputfile.open(filename, std::ios::binary);
      if (Inputfile.good()) {
	
	std::cout << "Binary reading "<<Global.size()*sizeof(double)
		  <<" bytes on "<< filename << " with offset "<< offset <<"\n";
	Inputfile.seekg(offset, std::ios::beg);
	Global.read_stream(Inputfile);
      }
      Inputfile.close();
    }

    int* global_site = new int[local_volume];
    std::vector< std::valarray<double>* > primaryField;//for accumulation
    std::valarray<double> local(format.size());
    
    for (int node=0; node < num_nodes; ++node) {
	primaryField.push_back(new std::valarray<double>(format.size()));
	for (int site = 0; site < local_volume; ++site)
	  {
	    global_site[site] = SiteIndex::instance()->
	      gsite(site);
	  }
	comm->sync();    
	//copy global index array on primary node
	comm->broadcast(local_volume, global_site[0], node);//1to1 enough but int not provided
	if(comm->primaryNode()) {	  
	  for (int site=0; site < format.Nvol(); ++site) 
	    {
	      for (int internal =0; internal < format.Nin(); ++internal) 
		{
		  for (int external =0; external < format.Nex(); ++external) 
		    {
		      local_index = format.index(internal, site, external);
		      //lime format: nb  external and site index are reverted from original order
		      global_index = internal +format.Nin()*(external+ format.Nex()*global_site[site]);
		      (*primaryField[node])[local_index] = Global[global_index];
		      
		    };
		};
	    };
	};
      }
	   //closing loop among processing nodes
	   comm->sync();
	 //copy one by one the local fields stored in the nodes into primary node
	 for (int node=0; node < num_nodes; ++node) 
       comm->send_1to1(local, *(primaryField[node]), local.size(),node, 0, node);	

    
     f = Field(local);

 delete[] global_site;
 
 return 0;


  }

   template <typename T>
     int ReadFromDisk(std::vector<Field>& f, const char* filename, int num_objects) {
     //this functions should be heavily modified
     T format(CommonPrms::instance()->Nvol());
     Field step(format.size());
     
     //informations about sizes should be provided by the file
      for (int field_num = 0; field_num < num_objects; ++field_num) {
	ReadFromDisk<T>(step, 
			filename, 
			sizeof(double)*field_num*format.size()*CommonPrms::instance()->NP());
	f.push_back(step);
      }
    }

}//end of namespace CCIO
