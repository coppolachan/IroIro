/*! 
 * @file fields_io_NERSC.cpp 
 *
 * @brief Declarations of MPI safeNERSC format reading routines for fields
 *
 * Time-stamp: <2013-06-03 16:00:50 neo>
 */

#include "fields_io.hpp"

#include <stdio.h>
#include <string.h>

namespace CCIO {
  Text_header::Text_header(){
  };

  void Text_header::add_element(std::string key, std::string value){
    tokens.insert(std::pair<std::string, std::string>(key, value));
  }

  void Text_header::get_value(std::string key, int& value){
    value = atoi(tokens[key].c_str());
  }

  void Text_header::get_value(std::string key, double& value){
    value = atof(tokens[key].c_str());
  }

  void Text_header::get_value(std::string key, std::string& value){
    value.assign(tokens[key].c_str());
  }

  void Text_header::print(){
    std::map<std::string, std::string>::iterator it;
    for (it=tokens.begin(); it!=tokens.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
  }

  //////////////////////////////////////////////////////////////////
  QCDheader* ReadNERSCheader(FILE *in, fpos_t *pos) {
#define MAX_LINE_LENGTH 1024
    Text_header* TH = new Text_header;
    char *key, *val;
    char *v;
    char line[MAX_LINE_LENGTH];
    int len;
    /* Begin reading, and check for "BEGIN_HEADER" token */
    v = fgets(line,MAX_LINE_LENGTH,in);
    if (strcmp(line,"BEGIN_HEADER\n")!=0){
      CCIO::cout << "NERSC type header not found!\nExiting...\n";
      exit(1);
    }
    
    while(1){
      v = fgets(line,MAX_LINE_LENGTH,in);
      if (strcmp(line,"END_HEADER\n")==0) {
	break;}

      /* Divide in tokens */
      key = strtok(line, " =");
      val = strtok(NULL, "\n");
      /* Store in the map for future search */
      TH->add_element(std::string(key), std::string(val+2));
    }
   
    return TH;
  }

  //////////////////////////////////////////////////////////////////
  void NERSCtoIROIRO_f(double* out, FILE* inputFile, int block, bool is3x2){
    //Allocate enough space
    float *uin = (float*) malloc(block*sizeof(float));
    size_t res = ReadStdBinary_Float(uin, inputFile, block);
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
  void NERSCtoIROIRO_d(double* out,FILE* inputFile, int block, bool is3x2){
    //Allocate enough space
    double *uin = (double*) malloc(block*sizeof(double));
    size_t res = ReadStdBinary(uin, inputFile, block);
#ifndef BIG_ENDIAN_TYPE
    byte_swap(uin, block);
#endif
    if (is3x2)
      reconstruct3x3<double>(out, uin, block);
    else
      memcpy(out, uin, sizeof(uin));
    free(uin);
  }
  
  ///////////////////////////////////////////////////////////////////////
  void ReadNERSCFormat(double* buffer,FILE *inputFile,
		       int block_size,
		       int tot_vol,int tot_in,int tot_ex,
		       QCDheader* QH){
    size_t res;
    std::string datatype;
    std::string floating_point; //32 or 64
    char * uin;
    bool is3x2 = false;
    QH->get_value("DATATYPE", datatype);
    QH->get_value("FLOATING_POINT", floating_point);
    
    int bl = 0;
    if (datatype.compare("4D_SU3_GAUGE")==0){
      bl= block_size/3*2;
      is3x2 = true;
    }
    else 
      if (datatype.compare("4D_SU3_GAUGE_3x3")==0)
	bl= block_size;
    
    if (bl == 0)
      Errors::IOErr(Errors::GenericError, "Corrupted header");
    
    if ((floating_point.compare("IEEE32")==0) ||(floating_point.compare("IEEE32BIG")==0) ){
      NERSCtoIROIRO_f(buffer, inputFile,  bl, is3x2);
    }
    else 
      if (floating_point.compare("IEEE64BIG")==0){
	NERSCtoIROIRO_d(buffer, inputFile, bl, is3x2);
      }
    
    // perform checks against header information
    /*
    QH->print();

    int d[4];
    QH->get_value("DIMENSION_1", d[0]);
    QH->get_value("DIMENSION_2", d[1]);
    QH->get_value("DIMENSION_3", d[2]);
    QH->get_value("DIMENSION_4", d[3]);
    std::ostringstream msg;
    if (d[0]!=CommonPrms::instance()->Lx()){
      msg << "Mismatch: XML input file  Lx = "<< CommonPrms::instance()->Lx()
	  << " - NERSC file  Lx = "<< d[0] << "\n";
      Errors::ParameterErr(msg);
    }
    if (d[1]!=CommonPrms::instance()->Ly()){
      msg << "Mismatch: XML input file  Ly = "<< CommonPrms::instance()->Ly()
	  << " - NERSC file  Ly = "<< d[1] << "\n";
      Errors::ParameterErr(msg);
    }
    if (d[2]!=CommonPrms::instance()->Lz()){
      msg << "Mismatch: XML input file  Lz = "<< CommonPrms::instance()->Lz()
	  << " - NERSC file  Lz = "<< d[2] << "\n";
      Errors::ParameterErr(msg);
    }
    if (d[3]!=CommonPrms::instance()->Lt()){
      msg << "Mismatch: XML input file  Lt = "<< CommonPrms::instance()->Lt()
	  << " - NERSC file  Lt = "<< d[3] << "\n";
      Errors::ParameterErr(msg);
    }

    */

  }


  ///////////////////////////////////////////////////////////////////////////
  QCDheader* NOheader(FILE *inputFile, fpos_t *pos){};
}
