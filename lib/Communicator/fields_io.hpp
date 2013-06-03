/*! 
 * @file fields_io.hpp 
 *
 * @brief Declarations of MPI safe read/write routines for fields
 *
 * Time-stamp: <2013-06-03 11:50:39 neo>
 */
#ifndef FIELDS_IO_HPP_
#define FIELDS_IO_HPP_

#include "include/field.h"
#include "include/messages_macros.hpp"
#include "include/errors.hpp"
#include "Geometry/siteIndex.hpp"
#include "Tools/byteswap.hpp"
#include <stdio.h>
#include <map>

#define CCIO_FILE_APPEND_MODE true
#define FORTRAN_CONTROL_WORDS 4  //number of fortran bytes for control

///////////////////////////////////////////////////////////////////
namespace CCIO {
  // definitions of storage format
  typedef int (*StorageFormat)(int, int, int, 
			       int, int, int);
  
  typedef void (*GenericReader)(double*,FILE *,int,
				int ,int ,int);

  //Abstract class for headers
  class QCDheader {
  public:
    virtual void get_value(std::string, int& ) = 0;
    virtual void get_value(std::string, double& ) = 0;
    virtual void get_value(std::string, std::string& ) = 0;
    
  };

  class Text_header: public QCDheader {
    std::map<std::string, char*> tokens;
  public:
    void add_element(std::string, char*);
    void get_value(std::string, int&);
    void get_value(std::string, double& );
    void get_value(std::string, std::string& );
    Text_header();

  };
  
  typedef QCDheader* (*GenericHeader)(FILE*, fpos_t *);

  typedef struct  {
    StorageFormat format;
    GenericReader reader;
    GenericHeader header;
    bool has_header;
  } ReaderFormat;
  


  //////////////////////////////////////////////////
  // ILDG binary storage format
  // same as NERSC (that eventually need reconstruction of third row)
  // same as MILC
  inline int ILDGBinFormat(int gsite,int in,int ex,
			   int tot_vol,int tot_in,int tot_ex){
    return in +tot_in*(ex+ tot_ex*gsite);
  }
  /////////////////////////////////////////////////
  inline int JLQCDLegacyFormat(int gsite,int in,int ex,
			       int tot_vol,int tot_in,int tot_ex){
    return in +tot_in*(gsite+ tot_vol*ex);
  }
  /////////////////////////////////////////////////  
  // reading function for corresponding format
  inline size_t ReadStdBinary(double* buffer,FILE *inputFile,
			      int block_size){
    //Read just as it is
    //order is (in, ex, sites)  
    size_t res = fread(buffer,sizeof(double),block_size,inputFile);
    return res;}
 
  inline size_t ReadStdBinary_Float(float* buffer,FILE *inputFile,
				    int block_size){
    //Read just as it is
    //order is (in, ex, sites)  
    size_t res = fread(buffer,sizeof(float),block_size, inputFile);
    return res;
} 
  
  inline void ReadILDGBFormat(double* buffer,FILE *inputFile,
			      int block_size,
			      int tot_vol,int tot_in,int tot_ex){
    size_t res = ReadStdBinary(buffer, inputFile, block_size);
  }
  
  void reconstruct3x3(double* out, float* in, int block_size);

  inline void ReadNERSCFormat(double* buffer,FILE *inputFile,
		       int block_size,
		       int tot_vol,int tot_in,int tot_ex){
    int bl = block_size/3*2;
    float* uin = (float *) malloc(bl*sizeof(float));
    size_t res = ReadStdBinary_Float(uin, inputFile, bl);//now assuming 3x2 format
    byte_swap(uin, bl);
    reconstruct3x3(buffer, uin, bl);
    // perform checks against header information
    free(uin);
  }

  QCDheader* NOheader(FILE *inputFile, fpos_t *pos);

  QCDheader* ReadNERSCheader(FILE *inputFile, fpos_t *pos);

  inline void ReadJLQCDLegacyFormat(double* buffer,FILE *inputFile,
				    int block_size,
				    int tot_vol,int tot_in,int tot_ex){
    //order is (in, sites, ex)
    //Read #tot_ex blocks
    int chunk = block_size/tot_ex;
    fpos_t pos;
    fgetpos(inputFile, &pos);
    for(int ext=0; ext<tot_ex; ++ext){
      //reads the ex blocks
      size_t res = fread((buffer+ext*chunk), sizeof(double), chunk, inputFile);
      fseek(inputFile, sizeof(double)*tot_vol*tot_in , SEEK_CUR);
    }
    fsetpos(inputFile, &pos);
    fseek(inputFile, sizeof(double)*chunk, SEEK_CUR);
  }


  ////////////////////////////////////////////////////////////////////

  /*!
   * @brief Saves a Field on the disk 
   * @param append_mode Used to write several fields in the same file (default = off)
   * @param T is the format specification
   */
  template <typename T>
  int SaveOnDisk(const Field& f, const char* filename, bool append_mode = false){
    Communicator* comm = Communicator::instance();
    CommonPrms* cmprms = CommonPrms::instance();
    T fmt(cmprms->Nvol());

    size_t local_idx, global_idx;
    int num_nodes, total_vol;
     
    vector_int block_x(cmprms->Nx());
    varray_double copy(fmt.Nin()*block_x.size()*fmt.Nex());
    varray_double local(fmt.Nin()*block_x.size()*fmt.Nex());

    // Open the output file, uses C style
    FILE * outFile;
    if(comm->primaryNode()) {
      if(append_mode) 	outFile = fopen (filename,"a");
      else          	outFile = fopen (filename,"w");

      std::cout << "Writing (binary mode) "<<f.size()*comm->size()*sizeof(double)
		<<" bytes on "<< filename << "... ";
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
		  for (int x_idx = 0; x_idx < block_x.size(); ++x_idx)
		    block_x[x_idx] =
		      SiteIndex::instance()->site(x_idx, y_slice, z_slice, t_slice);

		  copy = f[fmt.get_sub(block_x)];
	
		  // copy to the master node
		  comm->send_1to1(local, copy, local.size(), 0, node, node);
	
		  comm->sync();    
		  //Save to file sequentially
		  if (comm->primaryNode()){
		    for (int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for (int ex =0; ex < fmt.Nex(); ++ex) {	      
			for (int in =0; in < fmt.Nin(); ++in) {
			  //Breaks universality
			  local_idx = in +fmt.Nin()*(x_idx + block_x.size()*ex);
			  global_idx = in +fmt.Nin()*(ex+ fmt.Nex()*x_idx);//ILDG default
			  copy[global_idx] = local[local_idx];
			}
		      }
		    }
		    if (outFile!=NULL) 
		      fwrite((const char*)&copy[0], sizeof(double), 
			     local.size(), outFile);
		  }// end of if
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

  /*!
   * @brief Reads a Field from the disk 
   * @param storeFormat Function pointer to the storing format (ILDG / JLQCD legacy)
   * @param offset Starting point in file reading (byte offset)
   */
  template <typename T>
  int ReadFromDisk(Field& f,const char* filename,int offset, ReaderFormat readFormat){
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
    
    // Open the output file
    FILE * inFile;
    if(comm->primaryNode()){
      inFile = fopen(filename,"r");
      
      if (inFile == NULL)
	Errors::IOErr(Errors::FileNotFound, filename);

      std::cout << "Reading (binary mode) "<<f.size()*comm->size()*sizeof(double)
		<<" bytes from "<< filename << " with offset "<< offset <<"... ";
      fseek(inFile, offset, SEEK_SET);
    }
    
    if (readFormat.has_header){
      fpos_t pos;
      readFormat.header(inFile, &pos);
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
		      readFormat.reader((double*)&copy[0],inFile,
					copy.size(),total_vol,
					fmt.Nin(),fmt.Nex());
		      /*
		      if(storeFormat == ILDGBinFormat) 
			ReadILDGBFormat((double*)&copy[0],inFile,
					copy.size(),total_vol,
					fmt.Nin(),fmt.Nex());
		      if(storeFormat == JLQCDLegacyFormat) 
			ReadJLQCDLegacyFormat((double*)&copy[0],inFile,
					      copy.size(),total_vol,
					      fmt.Nin(),fmt.Nex());
		      */
		    }
		    for(int x_idx = 0; x_idx < block_x.size(); ++x_idx){
		      for(int ex =0; ex < fmt.Nex(); ++ex){	      
			for(int in =0; in < fmt.Nin(); ++in){
			  //Breaks universality
			  local_idx = in +fmt.Nin()*(x_idx +block_x.size()*ex);
			  /*
			  global_idx = storeFormat(x_idx,in,ex,block_x.size(),
						   fmt.Nin(),fmt.Nex());
			  */
			  global_idx = readFormat.format(x_idx,in,ex,block_x.size(),
						   fmt.Nin(),fmt.Nex());
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
    if(comm->primaryNode()){
      fclose(inFile);
      std::cout << "done\n";
    }
    return 0;
  }

  // With default reader
  template <typename T>
  int ReadFromDisk(Field& f,const char* filename, int offset = 0){
    ReaderFormat  DefaultStorage = {ILDGBinFormat,ReadILDGBFormat, NOheader, false};
    ReadFromDisk<T>(f, filename, offset, DefaultStorage);
  }
  
  template <typename T>
  int ReadFromDisk(std::vector<Field>& f, const char* filename, int num_objects) {
    //this functions should be heavily modified
    T fmt(CommonPrms::instance()->Nvol());
    Field step(fmt.size());
    
    //informations about sizes should be provided by the file
    for (int field_num = 0; field_num < num_objects; ++field_num) {
      ReadFromDisk<T>(step,filename,
		      sizeof(double)*field_num*fmt.size()*CommonPrms::instance()->NP());
      f.push_back(step);
    }
  }

}//end of namespace CCIO

#endif //FIELDS_IO_HPP
