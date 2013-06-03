/*! 
 * @file fields_io.cpp 
 *
 * @brief Declarations of MPI safe read/write routines for fields
 *
 * Time-stamp: <2013-06-03 11:52:49 neo>
 */

#include "fields_io.hpp"

#include <stdio.h>
#include <string.h>

namespace CCIO {
  Text_header::Text_header(){
  };

  void Text_header::add_element(std::string key, char* value){
    tokens.insert(std::pair<std::string, char*>(key, value));
  }

  void Text_header::get_value(std::string key, int& value){
    value = atoi(tokens[key]);
  }

  void Text_header::get_value(std::string key, double& value){
    value = atof(tokens[key]);
  }

  void Text_header::get_value(std::string key, std::string& value){
    value.assign(tokens[key]);
  }


  QCDheader* ReadNERSCheader(FILE *in, fpos_t *pos) {
#define MAX_LINE_LENGTH 1024
    char line[MAX_LINE_LENGTH];
    std::cout << "Header\n";
    /* Begin reading, and check for "BEGIN_HEADER" token */
    fgets(line,MAX_LINE_LENGTH,in);
    if (strcmp(line,"BEGIN_HEADER\n")!=0)
      exit(1);
    
    while(1){
      fgets(line,MAX_LINE_LENGTH,in);
      if (strcmp(line,"END_HEADER\n")==0) {
	std::cout << "Found end of header\n";
	break;}
    }
    

  }

  QCDheader* NOheader(FILE *inputFile, fpos_t *pos){};

  void reconstruct3x3(double* out, float* in,int block3x2){
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
