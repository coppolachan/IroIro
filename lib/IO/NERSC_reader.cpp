/*! 
 * @file NERSC_reader.cpp 
 *
 * @brief Definition of the NERSCReader class methods
 *
 * Time-stamp: <2013-06-04 17:58:33 neo>
 */
#include "NERSC_reader.hpp"
#include <string.h>
#include <stdlib.h>
#include "include/errors.hpp"
#include "Tools/byteswap.hpp"
#include "Measurements/GaugeM/staples.hpp"

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
    
    

   // check lattice size here
  }

  int NERSCReader::check(Field& U){
    CCIO::cout << "Checks against NERSC header\n";
    Staples stpl;
    double plaq = stpl.plaquette(GaugeField(U));
    CCIO::cout << "Plaquette: "<< plaq << " Stored value: "<< stored_plaquette;
    if (abs(plaq - stored_plaquette)< tolerance)
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
    else{
      for (int i =0; i < block; i++)
	out[i] = (double) uin[i];
    }
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
    else
      memcpy(out, uin, sizeof(uin));
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


