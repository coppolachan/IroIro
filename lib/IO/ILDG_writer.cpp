/*! 
 * @file ILDG_writer.cpp 
 *
 * @brief Definition of the ILDGWriter class methods
 *
 * Time-stamp: <2013-09-13 10:21:13 cossu>
 */

#include "ILDG_writer.hpp"
#include "include/messages_macros.hpp"
#include "Tools/byteswap.hpp"
#include <string.h>

namespace CCIO {
  void ILDGWriter::set_output(FILE *out){
    outputFile = out;
  }
  int ILDGWriter::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int ILDGWriter::write(double *buffer, unsigned int size){
    LimeRecordHeader *h;
    int status = EXIT_SUCCESS;
    if (writer_first_call){
      writer_first_call = 0;
 
      /* Write record header */
      CCIO::cout<< "Creating Header\n";
      h = limeCreateHeader(0, 1, "ildg-binary-data", size);
      
      CCIO::cout<< "Writing Header\n";
      status = limeWriteRecordHeader( h, LimeW );

      if( status < 0 ) { 
	CCIO::cerr<< "Header error\n";
	return EXIT_FAILURE;
      }
      
      limeDestroyHeader(h);
    }

    n_uint64_t bytes_to_write = size*sizeof(double);
    n_uint64_t bytes = bytes_to_write;
    CCIO::cout << "bytes "<<bytes<<"\n";
    status = limeWriteRecordData((char*)buffer, &bytes, LimeW);
    if( status != LIME_SUCCESS ) { 
      CCIO::cout << "LIME write error "<<status<<"\n";
      return EXIT_FAILURE;
    }
    
    num_blocks--;
    if (num_blocks == 0) {// last record
      limeWriterCloseRecord(LimeW);
    }
  }

  int ILDGWriter::header(){
    _Message(DEBUG_VERB_LEVEL, "ILDG Header...\n");
    n_uint64_t nbytes;
    int MB_flag=1, ME_flag=0;
    
    char message[] = "LIME test message";
    nbytes = strlen(message);

    LimeW = limeCreateWriter(outputFile);
    LimeHeader = limeCreateHeader( MB_flag, ME_flag, "lime-test-text", nbytes);
    limeWriteRecordHeader(LimeHeader, LimeW);
    limeDestroyHeader(LimeHeader);
    limeWriteRecordData(message, &nbytes, LimeW);
    limeWriterCloseRecord(LimeW);
    _Message(DEBUG_VERB_LEVEL, "ILDG Header: writing "<<nbytes<< " bytes record...\n");
  }

  ILDGWriter::ILDGWriter(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex),
    writer_first_call(1),
    num_blocks(tot_vol/CommonPrms::instance()->Nx()){};
}


