/*! 
 * @file fields_io.hpp 
 * @brief Definition of MPI safe input/output routines for fields
 *
 * These are the parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#ifndef FIELDS_IO_HPP_
#define FIELDS_IO_HPP_

#include "include/field.h"
#include "Main/Geometry/siteIndex.hpp"
#include "Tools/byteswap.hpp"

#ifdef IBM_BGQ
#include <stdio.h>
#endif

#define CCIO_FILE_APPEND_MODE true
#define FORTRAN_CONTROL_WORDS 4  //number of fortran bytes for control

namespace CCIO {
  
  inline int ILDGBinFormat(int gsite, 
			   int internal, 
			   int external, 
			   int tot_volume,
			   int tot_internal, 
			   int tot_external)  {
    return (internal +tot_internal*(external+ tot_external*gsite));
  }
  
  inline int JLQCDLegacyFormat(int gsite, 
			       int internal, 
			       int external, 
			       int tot_volume,
			       int tot_internal, 
			       int tot_external)  {
    return (internal +tot_internal*(gsite+ tot_volume*external));
  }
  
  typedef int (*StoringFormat)(int, int, int , int, int, int);

  /*!
   * @brief Saves a Field on the disk 
   *
   * @param append_mode Used to write several fields in the same file (default = off)
   */
  template <typename T>
  int SaveOnDisk(const Field& f, const char* filename, bool append_mode = false) 
  {
    T format(CommonPrms::instance()->Nvol());
    Field Global;
    Communicator* comm = Communicator::instance();
    size_t local_index, global_index;
    int num_nodes, total_volume;

    std::vector<int> global_site(format.Nvol());
    
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

      global_site = SiteIndex::instance()->get_gsite();
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
      #ifdef IBM_BGQ //hopefully temporary hack
      FILE * outFile;
      if(append_mode) 
	outFile = fopen (filename,"a");
      else
	outFile = fopen (filename,"w");

      if (outFile!=NULL) {
	std::cout << "Binary writing "<<f.size()*comm->size()*sizeof(double)
		  <<" bytes on "<< filename << " - BGQ hack\n";
	Global.write_buffer(outFile);
      }
      fclose(outFile);
      #else
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
      #endif
    }
    CCIO::cout << "Write ended succesfully\n";

    return 0;
  }
  
  /*!
   * @brief Saves several fields in the same file
   */
  template <typename T>
  int SaveOnDisk(const std::vector<Field>& f, const char* filename) {
    int result;
    for (int field_num = 0; field_num < f.size(); ++field_num) {
      result = SaveOnDisk<T>(f[field_num], filename, CCIO_FILE_APPEND_MODE);
    }
    return result;
  }

  
  /*!
   * @brief Reads a Field from the disk 
   *
   * @param storeFormat Function pointer to the storing format (ILDG / JLQCD legacy)
   * @param offset Starting point in file reading (byte offset)
   */

  template <typename T>
  int ReadFromDisk(Field& f, const char* filename, int offset = 0, StoringFormat storeFormat = ILDGBinFormat) {
    Communicator* comm = Communicator::instance();
    if(comm->primaryNode()) {std::cout << "Format initialization... " << std::endl;}
    T format(CommonPrms::instance()->Nvol());
    size_t local_index, global_index;
    int num_nodes, total_volume, local_volume;
    
    if(comm->primaryNode()) {std::cout << "Lattice constants initialization... " << std::endl;}
    total_volume = CommonPrms::instance()->Lvol();
    local_volume = CommonPrms::instance()->Nvol();
    num_nodes    = CommonPrms::instance()->NP();
    
    if(comm->primaryNode()) {std::cout << "Global field initialization... " << std::endl;}
    Field Global;


    //Read data from single file on disk
    if(comm->primaryNode()) {
      std::ifstream Inputfile;
      std::cout << "Global resize "  << format.size() << "*" << num_nodes << std::endl;
      Global.resize(format.size()*num_nodes);
      
      std::cout << "Opening file" << std::endl;
      Inputfile.open(filename, std::ios::binary);
      
      if (Inputfile.good()) {
	
	std::cout << "Binary reading "<<Global.size()*sizeof(double)
		  <<" bytes on "<< filename << " with offset "<< offset <<"... ";
	Inputfile.seekg(offset, std::ios::beg);
	if (storeFormat== ILDGBinFormat) {
	  Global.read_stream(Inputfile);
	} else {
	  if (storeFormat == JLQCDLegacyFormat) {
	    double* storage = new double[Global.size()];
	    std::cout << "JLQCDLegacyFormat ";
	    Inputfile.read((char*)storage, sizeof(double)*Global.size());
	    #ifdef BIG_ENDIAN_TYPE
	    //byte_swap_double(storage, Global.size());
	    #endif
	    Global = Field(std::valarray<double>(storage,Global.size())); 
	    delete[] storage;
	  }
	}
      } else {
	std::cout<< "Error in opening file [" << filename << "]\n";
	abort();
      }
      std::cout << "done\n";
      std::cout.flush();
      Inputfile.close();
    }
    
    std::vector<int> global_site(local_volume);
    //    std::vector< std::valarray<double>* > primaryField;//for accumulation
    std::valarray<double> local(format.size());
    std::valarray<double>* primaryField = new std::valarray<double>(format.size());   

    if(comm->primaryNode()) {std::cout << "Broadcasting on nodes... ";}
    for (int node=0; node < num_nodes; ++node) {
      //if(comm->primaryNode()) {std::cout << "node "<< node <<" starting loop \n";}
      //primaryField.push_back(new std::valarray<double>(format.size()));
      
      global_site = SiteIndex::instance()->get_gsite();
      comm->sync();    
      //copy global index array on primary node
      comm->broadcast(local_volume, global_site[0], node);//1to1 enough but int not provided
      //if(comm->primaryNode()) {std::cout << "node "<< node <<" broadcasting \n";}
      if(comm->primaryNode()) {	  
	for (int site=0; site < format.Nvol(); ++site) {
	  for (int internal =0; internal < format.Nin(); ++internal) {
	    for (int external =0; external < format.Nex(); ++external) {
	      local_index = format.index(internal, site, external);
	      global_index = storeFormat(global_site[site], internal, external, 
					 total_volume, format.Nin(), format.Nex());
	      //(*primaryField[node])[local_index] = Global[global_index];
	      (*primaryField)[local_index] = Global[global_index];
	      
	    }
	  }
	}
      } // end of index translation on primary Node
      comm->sync();
      comm->send_1to1(local, *primaryField, local.size(),node, 0, node);
      f = Field(local);
      //close node for loop
    }
    if(comm->primaryNode()) {std::cout << "done\n";}
    //closing loop among processing nodes - older version
    //comm->sync();
    //copy one by one the local fields stored in the nodes into primary node
    //for (int node=0; node < num_nodes; ++node) 
    //  comm->send_1to1(local, *(primaryField[node]), local.size(),node, 0, node);	
    
    //f = Field(local);
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


#endif //FIELDS_IO_HPP
