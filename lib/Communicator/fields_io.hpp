/*! 
 * @file fields_io.hpp 
 * @brief Definition of MPI safe read/write routines for fields
 *
 */

#ifndef FIELDS_IO_HPP_
#define FIELDS_IO_HPP_

#include "include/field.h"
#include "include/messages_macros.hpp"
#include "Main/Geometry/siteIndex.hpp"
#include "Tools/byteswap.hpp"
#include <stdio.h>

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


  
  inline void ReadILDGBFormat(double* buffer,
			     FILE *inputFile,
			     int block_size,
			     int tot_volume,
			     int tot_internal, 
			     int tot_external) {
    //Read just as it is
    //order is (internal, external, sites)  
    fread(buffer, sizeof(double), block_size, inputFile);
  }

  inline void ReadJLQCDLegacyFormat(double* buffer,
				    FILE *inputFile,
				    int block_size,
				    int tot_volume,
				    int tot_internal, 
				    int tot_external) {
    //order is (internal, sites, external)
    //Read #tot_external blocks
    int chunk = block_size/tot_external;
    fpos_t pos;
    fgetpos(inputFile, &pos);
    for (int ext=0; ext<tot_external; ++ext){
      //reads the external blocks
      fread((buffer+ext*chunk), sizeof(double), chunk, inputFile);
      fseek(inputFile, sizeof(double)*tot_volume*tot_internal , SEEK_CUR);
    }
    fsetpos(inputFile, &pos);
    fseek(inputFile, sizeof(double)*chunk, SEEK_CUR);
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
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    T format(cmprms->Nvol());

    size_t local_index, global_index;
    int num_nodes, total_volume;
     
    std::vector<int> block_x(cmprms->Nx());
    std::valarray<double> copy(format.Nin()*block_x.size()*format.Nex());
    std::valarray<double> local(format.Nin()*block_x.size()*format.Nex());

    // Open the output file
    FILE * outFile;
    if(comm->primaryNode()) {
      if(append_mode) 
	outFile = fopen (filename,"a");
      else
	outFile = fopen (filename,"w");

      std::cout << "Binary writing "<<f.size()*comm->size()*sizeof(double)
		<<" bytes on "<< filename << "... ";
    }
      
    //Loop among nodes (master node)   
    for (int node_t = 0; node_t < cmprms->NPEt(); ++node_t) {
      for (int t_slice = 0; t_slice < cmprms->Nt(); ++t_slice){
	
	for (int node_z = 0; node_z < cmprms->NPEz(); ++node_z) {
	  for (int z_slice = 0; z_slice <cmprms->Nz(); ++z_slice){
	    
	    for (int node_y = 0; node_y < cmprms->NPEy(); ++node_y) {
	      for (int y_slice = 0; y_slice < cmprms->Ny(); ++y_slice){
		
		for (int node_x = 0; node_x < cmprms->NPEx(); ++node_x) {
		  int node = comm->nodeid(node_x, node_y, node_z, node_t);
		  //Copy a block of dimension Nx into the master node and save

		  //Fill vector with indices
		  for (int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		    block_x[x_idx] = SiteIndex::instance()->site(x_idx, y_slice, z_slice, t_slice);
		  }
		  copy = f[format.get_sub(block_x)];
	
		  comm->send_1to1(local, copy, local.size(), 0, node, node);//copy to master node
	
		  comm->sync();    
		  //Save to file sequentially
		  if (comm->primaryNode()){
		    for (int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for (int external =0; external < format.Nex(); ++external) {	      
			for (int internal =0; internal < format.Nin(); ++internal) {
			  local_index = internal +format.Nin()*(x_idx + block_x.size()*external);//Breaks universality
			  global_index = internal +format.Nin()*(external+ format.Nex()*x_idx);//ILDG
			  copy[global_index] = local[local_index];
			}
		      }
		    }
		    if (outFile!=NULL) {
		      fwrite((const char*)&copy[0], sizeof(double), local.size(), outFile);
		    }
		  }
		  
		}
	      }
	      
	    }
	  }
	}
      }
    }
    comm->sync();    
    if(comm->primaryNode()) {
      fclose(outFile);
      std::cout << "write completed succesfully\n";
    }

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
    _Message(DEBUG_VERB_LEVEL, "Format initialization...\n");
    T format(CommonPrms::instance()->Nvol());
    size_t local_index, global_index;
    int num_nodes, total_volume, local_volume;
    
    _Message(DEBUG_VERB_LEVEL, "Lattice constants initialization...\n");
    total_volume = CommonPrms::instance()->Lvol();
    local_volume = CommonPrms::instance()->Nvol();
    num_nodes    = CommonPrms::instance()->NP();
 

       
    std::vector<int> block_x(CommonPrms::instance()->Nx());
    std::valarray<double> copy(format.Nin()*block_x.size()*format.Nex());
    std::valarray<double> local(format.Nin()*block_x.size()*format.Nex());
    
    // Open the output file
    FILE * inFile;
    if(comm->primaryNode()) {
	inFile = fopen (filename,"r");
	
	std::cout << "Binary reading "<<f.size()*comm->size()*sizeof(double)
		  <<" bytes from "<< filename << " with offset "<< offset <<"... ";
	fseek(inFile, offset, SEEK_SET);
    }
    

    //Loop among nodes (master node)   
    for (int node_t = 0; node_t < CommonPrms::instance()->NPEt(); ++node_t) {
      for (int t_slice = 0; t_slice < CommonPrms::instance()->Nt(); ++t_slice){
	
	for (int node_z = 0; node_z < CommonPrms::instance()->NPEz(); ++node_z) {
	  for (int z_slice = 0; z_slice < CommonPrms::instance()->Nz(); ++z_slice){
	    
	    for (int node_y = 0; node_y < CommonPrms::instance()->NPEy(); ++node_y) {
	      for (int y_slice = 0; y_slice < CommonPrms::instance()->Ny(); ++y_slice){
		
		for (int node_x = 0; node_x < CommonPrms::instance()->NPEx(); ++node_x) {
		  int node = comm->nodeid(node_x, node_y, node_z, node_t);
		  //Copy a block of dimension Nx into the master node and save

		  //Read from file sequentially
		  if (comm->primaryNode()){
		    if (inFile!=NULL) {
		      if (storeFormat == JLQCDLegacyFormat) {
			ReadJLQCDLegacyFormat((double*)&copy[0], inFile, copy.size(), total_volume, format.Nin(), format.Nex());
		      }
		      if (storeFormat == ILDGBinFormat) {
			ReadILDGBFormat((double*)&copy[0], inFile, copy.size(), total_volume, format.Nin(), format.Nex());
		      }
		      
		    }

		    for (int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for (int external =0; external < format.Nex(); ++external) {	      
			for (int internal =0; internal < format.Nin(); ++internal) {
			  local_index = internal +format.Nin()*(x_idx + block_x.size()*external);//Breaks universality
			  global_index = storeFormat(x_idx, internal, external, 
						     block_x.size(), format.Nin(), format.Nex());
			  local[local_index] = copy[global_index];
			}
		      }
		    }
		  }//end primaryNode
		  comm->sync();    
		  comm->send_1to1(copy, local, local.size(), node, 0, node);//copy to master node
		  comm->sync();    
		  if(comm->nodeid() == node){
		    //Fill vector with indices
		    for (int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      block_x[x_idx] = SiteIndex::instance()->site(x_idx, y_slice, z_slice, t_slice);
		    }
		    
		    f.set(format.get_sub(block_x),copy);
		    
		  }


		  
		}
	      }
	      
	    }
	  }
	}
      }
    }
    comm->sync();    
    if(comm->primaryNode()) {
      fclose(inFile);
      std::cout << "done\n";
    }

    /*
    //Read data from single file on disk
    Field Global;
    if(comm->primaryNode()) {
      std::ifstream Inputfile;
      std::cout << "Global resize "  << format.size() << "*" << num_nodes << std::endl;
      Global.resize(format.size()*num_nodes);
      
      std::cout << "Opening file" << std::endl;
      Inputfile.open(filename, std::ios::binary);
      
      if (Inputfile.good()) {
	
	std::cout << "OLD - Binary reading "<<Global.size()*sizeof(double)
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
    */



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
