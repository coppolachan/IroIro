/*! 
 * @file NERSC_reader.cpp 
 *
 * @brief Definition of the NERSCReader class methods
 *
 * Time-stamp: <2014-07-15 14:28:42 neo>
 */
#include "NERSC_reader.hpp"
#include <string.h>
#include <stdlib.h>
#include "include/errors.hpp"
#include "Tools/byteswap.hpp"
#include "Measurements/GaugeM/staples.hpp"

#include <iostream>     // std::cout, std::hex, std::endl
#include <iomanip>      // std::setiosflags, std::resetiosflags

#define MAX_LINE_LENGTH 1024

namespace CCIO {
  void NERSCReader::set_sources(FILE *src){
    inputFile = src;
  }
  int NERSCReader::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int NERSCReader::read(double *buffer, unsigned int size){
    //Read just as it is
    //order is (in, ex, sites)  
    
    int bl = 0;
    if (is3x2)
      bl= size/3*2;
    else 
      bl= size;
    
    if (bl == 0)
      Errors::IOErr(Errors::GenericError, "Corrupted header");
    
    if (bits == 32)
      NERSCtoIROIRO_f(buffer, bl);
    else 
      NERSCtoIROIRO_d(buffer, bl);
  }

  int NERSCReader::header(){
    char *key, *val;
    char *v;
    char line[MAX_LINE_LENGTH];
    int len;
    /* Begin reading, and check for "BEGIN_HEADER" token */
    v = fgets(line,MAX_LINE_LENGTH,inputFile);
    if (strcmp(line,"BEGIN_HEADER\n")!=0){
      Errors::IOErr(Errors::GenericError, 
		    "NERSC type header not found!\nExiting...\n");
    }
    
    while(1){
      v = fgets(line,MAX_LINE_LENGTH,inputFile);
      if (strcmp(line,"END_HEADER\n")==0) {
	break;}

      /* Divide in tokens */
      key = strtok(line, " =");
      val = strtok(NULL, "\n");
      /* Store in the map for future search */
      TH.add_element(std::string(key), std::string(val+2));
    }

    std::string datatype;
    std::string floating_point; //32 or 64
    TH.get_value("DATATYPE", datatype);
    TH.get_value("FLOATING_POINT", floating_point);
    TH.get_value("PLAQUETTE", stored_plaquette);
    TH.get_value("LINK_TRACE", stored_link);

    


    if (datatype.compare("4D_SU3_GAUGE")==0)
      is3x2 = true;
    else
      if (datatype.compare("4D_SU3_GAUGE_3x3")==0)
	is3x2 = false;

    if ((floating_point.compare("IEEE32")==0) ||(floating_point.compare("IEEE32BIG")==0) ){
      bits = 32;
      tolerance = 1e-8;
    }
    else 
      if (floating_point.compare("IEEE64BIG")==0){
	bits = 64;
	tolerance = 1e-15;
      }
    
    CCIO::cout << "\n";
    CCIO::cout << "3x2 format : " << is3x2 << "\n";
    CCIO::cout << "Floating point format : " << bits << "-bit\n";
    
    // check lattice size here!!!!!!!
    int nersc_dim[4];
    TH.get_value("DIMENSION_1", nersc_dim[0]);
    TH.get_value("DIMENSION_2", nersc_dim[1]);
    TH.get_value("DIMENSION_3", nersc_dim[2]);
    TH.get_value("DIMENSION_4", nersc_dim[3]);

    if (nersc_dim[0] != CommonPrms::instance()->Lx()) 
      Errors::IOErr(Errors::GenericError, "NERSC file dimension X do not match the value in the XML file."); 

    if (nersc_dim[1] != CommonPrms::instance()->Ly()) 
      Errors::IOErr(Errors::GenericError, "NERSC file dimension Y do not match the value in the XML file."); 

    if (nersc_dim[2] != CommonPrms::instance()->Lz()) 
      Errors::IOErr(Errors::GenericError, "NERSC file dimension Z do not match the value in the XML file."); 

    if (nersc_dim[3] != CommonPrms::instance()->Lt()) 
      Errors::IOErr(Errors::GenericError, "NERSC file dimension T do not match the value in the XML file."); 

    CCIO::cout << "Dimensions check: OK\n";

    // checksum
    //save location in the file
    long unsigned int data_base_position = ftell(inputFile);

    int matrix_size;
    if(is3x2)
      matrix_size = 12;
    else
      matrix_size = 18;


    if (bits == 32){
      uint32_t stored_chksum = 0, chksum = 0;
      TH.get_value("CHECKSUM", stored_chksum); 
      
      //expected data size
      long unsigned int total_size = total_volume*total_external*matrix_size;
      long unsigned int chunk_size = total_size/nersc_dim[0];

      float *chk_buf = new float[chunk_size]; 
      for (int chunk = 0; chunk< nersc_dim[0]; chunk++){
	size_t res = fread(chk_buf, sizeof(float), chunk_size, inputFile);
	uint32_t* chk_ptr = (uint32_t*)chk_buf;
	for(unsigned int i=0; i < chunk_size*sizeof(float)/sizeof(uint32_t); i++)
	  chksum += chk_ptr[i];

      } 
      printf("Checksum: %x  Stored value: %x\n", chksum, stored_chksum );
      delete[] chk_buf;
    }
    if (bits == 64){
      uint32_t stored_chksum = 0, chksum = 0;
      TH.get_value("CHECKSUM", stored_chksum); 
      
      //expected data size
      long unsigned int total_size = total_volume*total_external*matrix_size;
      long unsigned int chunk_size = total_size/nersc_dim[0];

      double *chk_buf = new double[chunk_size]; 
      for (int chunk = 0; chunk< nersc_dim[0]; chunk++){
	size_t res = fread(chk_buf, sizeof(double), chunk_size, inputFile);
	uint32_t* chk_ptr = (uint32_t*)chk_buf;
	for(unsigned int i=0; i < chunk_size*sizeof(double)/sizeof(uint32_t); i++)
	  chksum += chk_ptr[i];

      }  
      printf("Checksum: %x  Stored value: %x\n", chksum, stored_chksum );
      delete[] chk_buf;
    }
    

    //reset location
    fseek(inputFile, data_base_position, SEEK_SET);
  }

  int NERSCReader::check(Field& U){
    CCIO::cout << "Checks against NERSC header\n";
    Staples stpl;
    double plaq = stpl.plaquette(GaugeField(U));
    CCIO::cout << "Plaquette: "<< plaq << " Stored value: "<< stored_plaquette << " Difference: "<< fabs(plaq - stored_plaquette) << " Tolerance: "<< tolerance << "\n";
    
    if (fabs(plaq - stored_plaquette)< tolerance)
      CCIO::cout << " (OK)\n";
    else
     Errors::IOErr(Errors::GenericError, "Something was wrong in the reading process."); 
  }

  NERSCReader::NERSCReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};

  //////////////////////////////////////////////////////////////////
  void NERSCReader::NERSCtoIROIRO_f(double* out,int block){
    //Allocate enough space
    float *uin = (float*) malloc(block*sizeof(float));
    size_t res = fread(uin,sizeof(float),block, inputFile);
#ifndef BIG_ENDIAN_TYPE
    byte_swap(uin, block);
#endif
    if (is3x2)
      reconstruct3x3<float>(out, uin, block);
    else
      memcpy(out, uin, block*sizeof(float));
    free(uin);
  }
  
  ////////////////////////////////////////////////////////////////////
  void NERSCReader::NERSCtoIROIRO_d(double* out, int block){
    //Allocate enough space
    double *uin = (double*) malloc(block*sizeof(double));
    size_t res = fread(uin,sizeof(double),block,inputFile);
    
#ifndef BIG_ENDIAN_TYPE
    byte_swap(uin, block);
#endif
    
    if (is3x2)
      reconstruct3x3<double>(out, uin, block);
    else{
      //for (int i =0; i < block; i++)
      //	out[i] = (double) uin[i];
      memcpy(out, uin, block*sizeof(double));
    }
    free(uin);
  }

template <typename FLOAT>
void NERSCReader::reconstruct3x3(double* out, FLOAT* in,int block3x2){
  int matrix_block3x2 = 12; //3x2x2
  int matrix_block3x3 = 18; //3x3x2
  double* _p;
  _p = out;
  for (int i =0; i < block3x2/matrix_block3x2; i++){
    //Copy first two columns
    for (int s =0; s < matrix_block3x2; s++){
      _p[s] = (double) in[i*matrix_block3x2+s];
    }
    _p[12] = _p[ 2]*_p[10] - _p[ 4]*_p[ 8] - _p[ 3]*_p[11] + _p[ 5]*_p[ 9];
    _p[13] = _p[ 4]*_p[ 9] - _p[ 2]*_p[11] + _p[ 5]*_p[ 8] - _p[ 3]*_p[10];
    _p[14] = _p[ 4]*_p[ 6] - _p[ 0]*_p[10] - _p[ 5]*_p[ 7] + _p[ 1]*_p[11];
    _p[15] = _p[ 0]*_p[11] - _p[ 4]*_p[ 7] + _p[ 1]*_p[10] - _p[ 5]*_p[ 6];
    _p[16] = _p[ 0]*_p[ 8] - _p[ 2]*_p[ 6] - _p[ 1]*_p[ 9] + _p[ 3]*_p[ 7];
    _p[17] = _p[ 2]*_p[ 7] - _p[ 0]*_p[ 9] + _p[ 3]*_p[ 6] - _p[ 1]*_p[ 8];	
    
    _p += matrix_block3x3;	
  }
}


}


