/*! 
 * @file fields_io.hpp 
 *
 * @brief Declarations of MPI safe read/write routines for fields
 *
 * Time-stamp: <2014-07-18 16:52:02 noaki>
 */
#ifndef FIELDS_IO_HPP_
#define FIELDS_IO_HPP_

#include "include/config.h"
#include "include/field.h"
#include "include/messages_macros.hpp"
#include "include/errors.hpp"
#include "Geometry/siteIndex.hpp"
#include "Geometry/siteIndex3d.hpp"
#include <stdio.h>
#include <string.h>

#include "generic_header.hpp"
#include "readers.hpp"
#include "writers.hpp"


#define CCIO_FILE_APPEND_MODE true
#define FORTRAN_CONTROL_WORDS 4  //number of fortran bytes for control

namespace CCIO {
  /*!
   * @brief Factory for Readers
   * @param ID string containing the ID of the input mode
   */
  GeneralReader* getReader(const std::string ID, int, int, int);

  /*!
   * @brief Factory for Readers
   * @param ID string containing the ID of the input mode
   */
  GeneralWriter* getWriter(const std::string ID, int, int, int, const Field&);
  GeneralWriter* getWriter(const std::string ID, int, int, int);

  /*!
   * @brief Saves a Field on the disk 
   * @param append_mode Used to write several fields in the same file (default = off)
   * @param T is the format specification
   */
  template <typename T>
  int SaveOnDisk(const Field& f, const char* filename,
		 bool append_mode = false,
		 const std::string writerID = "Binary"){
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    T fmt(cmprms->Nvol());
    GeneralWriter *writer;

    size_t local_idx, global_idx;
    int num_nodes, local_vol, total_vol;
    
    total_vol = cmprms->Lvol();
    local_vol = cmprms->Nvol();
    num_nodes = cmprms->NP();
    
    vector_int block_x(cmprms->Nx());
    varray_double copy(fmt.Nin()*block_x.size()*fmt.Nex());
    varray_double local(fmt.Nin()*block_x.size()*fmt.Nex());

    writer = getWriter(writerID, total_vol, fmt.Nin(), fmt.Nex(), f); 

    // Open the output file, uses C style
    FILE * outFile;
    if(comm->primaryNode()) {
      if(append_mode) 	outFile = fopen (filename,"a");
      else          	outFile = fopen (filename,"w");

      if(outFile == NULL) Errors::IOErr(Errors::FileNotFound, filename);

      std::cout << "Writing file of type ["<<writerID<<"] - size "
		<< f.size()*comm->size()*sizeof(double)
		<< " bytes on "<< filename << "... ";
      
      writer->set_output(outFile);
      writer->header();      
    }

    //Loop among nodes (master node)   
    for(int node_t = 0; node_t < cmprms->NPEt(); ++node_t) {
      for(int t_slice = 0; t_slice < cmprms->Nt(); ++t_slice){
	
	for(int node_z = 0; node_z < cmprms->NPEz(); ++node_z) {
	  for(int z_slice = 0; z_slice <cmprms->Nz(); ++z_slice){
	    
	    for(int node_y = 0; node_y < cmprms->NPEy(); ++node_y) {
	      for(int y_slice = 0; y_slice < cmprms->Ny(); ++y_slice){
		
		for(int node_x = 0; node_x < cmprms->NPEx(); ++node_x) {
		  int node = comm->nodeid(node_x, node_y, node_z, node_t);
		  //Copy a block of dimension Nx into the master node and save

		  //Fill vector with indices
		  for(int x_idx = 0; x_idx < block_x.size(); ++x_idx)
		    block_x[x_idx] =
		      SiteIndex::instance()->site(x_idx, y_slice, z_slice, t_slice);

		  copy = f[fmt.get_sub(block_x)];
	
		  // copy to the master node
		  comm->send_1to1(local,copy,local.size(),0,node,node);
	
		  comm->sync();    
		  //Save to file sequentially
		  if(comm->primaryNode()){
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for(int ex =0; ex < fmt.Nex(); ++ex) {	      
			for(int in =0; in < fmt.Nin(); ++in) {
			  //Breaks universality
			  local_idx = in +fmt.Nin()*(x_idx + block_x.size()*ex);
			  global_idx = writer->format(x_idx, in, ex);
			  copy[global_idx] = local[local_idx];
			}
		      }
		    }
		    if(outFile!=NULL) writer->write((double*)&copy[0], copy.size());
		  }// end of primaryNode
		}// end of node_x
	      }
	    }// end of node_y
	  }
	}// end of node_z
      }
    }// end of node_t

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
    for (int field_num = 0; field_num < f.size(); ++field_num) 
      result = SaveOnDisk<T>(f[field_num], filename, CCIO_FILE_APPEND_MODE);
    return result;
  }

  template <typename T>
  int SaveOnDisk3D(const Field& f, const char* filename,
		   bool append_mode = false,
		   const std::string writerID = "Binary"){
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    assert(cmprms->NPEt()==1); // only works with NPEt=1

    T fmt(cmprms->Nvol());
    GeneralWriter *writer;

    size_t local_idx, global_idx;
    int num_nodes, local_vol, total_vol;
    
    total_vol = cmprms->Lvol()/cmprms->Lt();
    local_vol = cmprms->Nvol()/cmprms->Nt();
    num_nodes = cmprms->NP();
    
    SiteIndex3d idx3d;

    vector_int block_x(cmprms->Nx());
    varray_double copy(fmt.Nin()*block_x.size()*fmt.Nex());
    varray_double local(fmt.Nin()*block_x.size()*fmt.Nex());

    writer = getWriter(writerID,total_vol,fmt.Nin(),fmt.Nex()); 

    // Open the output file, uses C style
    FILE * outFile;
    if(comm->primaryNode()) {
      if(append_mode) 	outFile = fopen (filename,"a");
      else          	outFile = fopen (filename,"w");

      if(outFile == NULL) Errors::IOErr(Errors::FileNotFound, filename);

      std::cout << "Writing file of type ["<<writerID<<"] - size "
		<< f.size()*comm->size()*sizeof(double)
		<< " bytes on "<< filename << "... ";
      
      writer->set_output(outFile);
      writer->header();      
    }

    //Loop among nodes (master node)   
    int node_t = 0;
    for(int node_z = 0; node_z < cmprms->NPEz(); ++node_z) {
      for(int z_slice = 0; z_slice <cmprms->Nz(); ++z_slice){
	    
	for(int node_y = 0; node_y < cmprms->NPEy(); ++node_y) {
	  for(int y_slice = 0; y_slice < cmprms->Ny(); ++y_slice){
		
	    for(int node_x = 0; node_x < cmprms->NPEx(); ++node_x) {
	      int node = comm->nodeid(node_x,node_y,node_z,node_t);
	      //Copy a block of dimension Nx into the master node and save

	      //Fill vector with indices
	      for(int x_idx = 0; x_idx < block_x.size(); ++x_idx)
		block_x[x_idx] = idx3d.site(x_idx, y_slice, z_slice);

	      copy = f[fmt.get_sub(block_x)];
	
	      // copy to the master node
	      comm->send_1to1(local,copy,local.size(),0,node,node);
	
	      comm->sync();    
	      //Save to file sequentially
	      if(comm->primaryNode()){
		for(int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		  for(int ex =0; ex < fmt.Nex(); ++ex) {	      
		    for(int in =0; in < fmt.Nin(); ++in) {
		      //Breaks universality
		      local_idx = in +fmt.Nin()*(x_idx + block_x.size()*ex);
		      global_idx = writer->format(x_idx, in, ex);
		      copy[global_idx] = local[local_idx];
		    }
		  }
		}
		if(outFile!=NULL) writer->write((double*)&copy[0], copy.size());
	      }// end of primaryNode
	    }// end of node_x
	  }
	}// end of node_y
      }
    }// end of node_z
    
    comm->sync();    
    if(comm->primaryNode()) {
      fclose(outFile);
      std::cout << "write completed succesfully\n";
    }
    return 0;
  }

  /////////////////////////////////////////////////////////////////////////////////////
  /*!
   * @brief Reads a Field from the disk 
   * @param storeFormat Function pointer to the storing format (ILDG / JLQCD legacy)
   * @param offset Starting point in file reading (byte offset)
   */
  template <typename T>
  int ReadFromDisk(Field& f,
		   const char* filename,
		   const uint64_t offset = 0, 
		   const std::string readerID = "Binary",
		   const bool do_check = false){
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    _Message(DEBUG_VERB_LEVEL, "Format initialization...\n");
    
    T fmt(cmprms->Nvol());
    
    size_t local_idx, global_idx;
    int num_nodes, total_vol, local_vol;
    
    _Message(DEBUG_VERB_LEVEL, "Lattice constants initialization...\n");
    
    total_vol = cmprms->Lvol();
    local_vol = cmprms->Nvol();
    num_nodes    = cmprms->NP();
       

    vector_int block_x(cmprms->Nx());
    varray_double copy(fmt.Nin()*block_x.size()*fmt.Nex());
    varray_double local(fmt.Nin()*block_x.size()*fmt.Nex());
    
    // Get the reader class
    GeneralReader* reader = getReader(readerID, total_vol, fmt.Nin(), fmt.Nex());

    // Open the output file
    FILE * inFile;
    if(comm->primaryNode()){
      inFile = fopen(filename,"r");
      
      if(inFile == NULL) Errors::IOErr(Errors::FileNotFound, filename);

      //check file size

      std::cout << "Reading file of type ["<<readerID<<"] - size "<<f.size()*comm->size()*sizeof(double)
		<<" bytes from "<< filename << " with offset "<< offset <<"... " << std::flush;
      fseek(inFile, offset, SEEK_SET);
    
      reader->set_sources(inFile);
      reader->header();

    }

    //Loop among nodes (master node)   
    for(int node_t = 0; node_t < cmprms->NPEt(); ++node_t){
      for(int t_slice = 0; t_slice < cmprms->Nt(); ++t_slice){
	
	for(int node_z = 0; node_z < cmprms->NPEz(); ++node_z){
	  for(int z_slice = 0; z_slice < cmprms->Nz(); ++z_slice){
	    
	    for(int node_y = 0; node_y < cmprms->NPEy(); ++node_y){
	      for(int y_slice = 0; y_slice < cmprms->Ny(); ++y_slice){
		
		for(int node_x = 0; node_x < cmprms->NPEx(); ++node_x){
		  int node = comm->nodeid(node_x, node_y, node_z, node_t);
		  //Copy a block of dimension Nx into the master node and save

		  //Read from file sequentially
		  if(comm->primaryNode()){
		    if(inFile!=NULL){
		      reader->read((double*)&copy[0],copy.size());
		    }
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for(int ex =0; ex < fmt.Nex(); ++ex){	      
			for(int in =0; in < fmt.Nin(); ++in){
			  //Breaks universality (check this)
			  local_idx = in +fmt.Nin()*(x_idx +block_x.size()*ex);
			  global_idx = reader->format(x_idx, in, ex);
			  local[local_idx] = copy[global_idx];
			}
		      }
		    }
		  }//end primaryNode
		  comm->sync();    
		  //copy to master node
		  comm->send_1to1(copy,local,local.size(),node,0,node);
		  comm->sync();    
		  if(comm->nodeid() == node){
		    //Fill vector with indices
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx)
		      block_x[x_idx] = SiteIndex::instance()->site(x_idx, y_slice,
								   z_slice, t_slice);
		    f.set(fmt.get_sub(block_x),copy);
		  }
		}
	      }
	    }// end of node_y
	  }
	}// end of node_z
      }
    }// end of node_t
    comm->sync();  

    if(do_check){
      if(reader->check(f) != CHECK_PASS) 
	Errors::IOErr(Errors::GenericError, 
		      "Failed some basic tests on the gauge configuration");
      CCIO::cout << "Tests passed!\n";
    }

    if(comm->primaryNode()){
      fclose(inFile);
      std::cout << "done\n";
    }
    return 0;
  }

  template <typename T>
  int ReadFromDisk(std::vector<Field>& f, const char* filename, int num_objects) {
    //this functions should be heavily modified
    T fmt(CommonPrms::instance()->Nvol());
    Field step(fmt.size());
    
    //informations about sizes should be provided by the file
    for (int field_num = 0; field_num < num_objects; ++field_num) {
      uint64_t offset = sizeof(double)*field_num*fmt.size()*CommonPrms::instance()->NP();
      ReadFromDisk<T>(step,filename,offset);
      f.push_back(step);
    }
  }



  /////////////////////////////////////////////////////////////////////////////////////
  /*!
   * @brief Reads a 3D Field from the disk. Need a specialized function
   * @param storeFormat Function pointer to the storing format (ILDG / JLQCD legacy)
   * @param offset Starting point in file reading (byte offset)
   */
  template <typename T>
  int ReadFromDisk3D(Field& f,
		     const char* filename,
		     const int block3D = 0,
		     const uint64_t offset = 0, 
		     const std::string readerID = "Binary",
		     const bool do_check = true){
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    _Message(DEBUG_VERB_LEVEL, "Format initialization...\n");
    
    T fmt(cmprms->Nvol());
    
    size_t local_idx, global_idx;
    int num_nodes, total_vol, local_vol;
    
    _Message(DEBUG_VERB_LEVEL, "Lattice constants initialization...\n");
    
    total_vol = cmprms->Lvol()/cmprms->Lt();
    local_vol = cmprms->Nvol()/cmprms->Nt();
    num_nodes    = cmprms->NP();
    
    
    vector_int block_x(cmprms->Nx());
    varray_double copy(fmt.Nin()*block_x.size()*fmt.Nex());
    varray_double local(fmt.Nin()*block_x.size()*fmt.Nex());
    
    // Get the reader class
    GeneralReader* reader = getReader(readerID, total_vol, fmt.Nin(), fmt.Nex());
    
    // Open the output file
    FILE * inFile;
    if(comm->primaryNode()){
      inFile = fopen(filename,"r");
      
      if(inFile == NULL) Errors::IOErr(Errors::FileNotFound, filename);
      
      //check file size
      int bytes_block3D = f.size()*comm->size()/CommonPrms::instance()->NPEt()*sizeof(double);

      std::cout << "Reading 3D file of type ["<<readerID<<"] - size "<<bytes_block3D
		<<" bytes from "<< filename << " with offset "<< offset <<"... " << std::flush;
      
      // Jump to the block3D-th record
      uint64_t block_offset = offset + block3D*bytes_block3D;

      fseek(inFile, block_offset, SEEK_SET);
      
      reader->set_sources(inFile);
      reader->header();
      
    }
    
    //Loop among nodes (master node) 
    int node_t = 0;
    //    for(int node_t = 0; node_t < cmprms->NPEt(); ++node_t){
    //for(int t_slice = 0; t_slice < blocks3D/(cmprms->NPEt()); ++t_slice){ // read as 4d field with T dim = 3d_blocks
    int t_slice = 0; // store in the first slice

	  
	  for(int node_z = 0; node_z < cmprms->NPEz(); ++node_z){
	    for(int z_slice = 0; z_slice < cmprms->Nz(); ++z_slice){
	      
	      for(int node_y = 0; node_y < cmprms->NPEy(); ++node_y){
	      for(int y_slice = 0; y_slice < cmprms->Ny(); ++y_slice){
		
		for(int node_x = 0; node_x < cmprms->NPEx(); ++node_x){
		  int node = comm->nodeid(node_x, node_y, node_z, node_t);
		  //Copy a block of dimension Nx into the master node and save
		  
		  //Read from file sequentially
		  if(comm->primaryNode()){
		    if(inFile!=NULL){
		      reader->read((double*)&copy[0],copy.size());
		    }
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for(int ex =0; ex < fmt.Nex(); ++ex){	      
			for(int in =0; in < fmt.Nin(); ++in){
			  //Breaks universality (check this)
			  local_idx = in +fmt.Nin()*(x_idx +block_x.size()*ex);
			  global_idx = reader->format(x_idx, in, ex);
			  local[local_idx] = copy[global_idx];
			}
		      }
		    }
		  }//end primaryNode
		  comm->sync();    
		  //copy to master node
		  comm->send_1to1(copy,local,local.size(),node,0,node);
		  comm->sync();    
		  if(comm->nodeid() == node){
		    //Fill vector with indices
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx)
		      block_x[x_idx] = SiteIndex::instance()->site(x_idx, y_slice,
								   z_slice, t_slice);
		    f.set(fmt.get_sub(block_x),copy);
		  }
		}
	      }
	    }// end of node_y
	  }
	}// end of node_z
    //  }
  //    }// end of node_t
    comm->sync();  

    if(do_check){
      if(reader->check(f) != CHECK_PASS) 
	Errors::IOErr(Errors::GenericError, 
		      "Failed some basic tests on the gauge configuration");
      CCIO::cout << "Tests passed!\n";
    }

    if(comm->primaryNode()){
      fclose(inFile);
      std::cout << "done\n";
    }
    return 0;
  }

}//end of namespace CCIO

#endif //FIELDS_IO_HPP
