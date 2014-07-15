/*! 
 * @file NERSC_writer.cpp 
 *
 * @brief Definition of the NERSCWriter class methods
 *
 * Time-stamp: <2014-07-15 14:33:10 neo>
 */

#include "NERSC_writer.hpp"
#include "include/messages_macros.hpp"
#include "Tools/byteswap.hpp"
#include "include/errors.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include <string.h>
#include <stdint.h>

namespace CCIO {
  unsigned int NERSCWriter::checksum(const GaugeField &U, int matrix_block ){
    size_t mat_size = sizeof(double)*matrix_block; //matrix_block = 12,18
    uint32_t checksum = 0;   // checksum
    //Allocate space for array
    //char *chk_buf = new char[mat_size];
    double *chk_buf = new double[matrix_block];

    //Get pointer to the block of data
    for (int s = 0; s < CommonPrms::instance()->Nvol(); s++){
      for (int dir = 0; dir < NDIM_; dir++){
	double* base_ptr = const_cast<GaugeField&>(U).data.getaddr(U.format.index(0,s,dir));
	memcpy(chk_buf, base_ptr, mat_size);
	
	//Eventually flip bits
#ifndef BIG_ENDIAN_TYPE
	byte_swap_double(chk_buf, mat_size);
#endif
	//Compute checksum
	uint32_t* chk_ptr = (uint32_t*)chk_buf;
	for(unsigned int i=0; i < mat_size/sizeof(uint32_t); ++i)
	  checksum += chk_ptr[i];

      }
    }
    //Delete pointer
    delete[] chk_buf;

    //Global sum ints
    return Communicator::instance()->reduce_sum(checksum);
  }

  void NERSCWriter::set_output(FILE *src){
    outputFile = src;
  }
  int NERSCWriter::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int NERSCWriter::header() {
    CCIO::cout << "Writing NERSC header\n";
    // write in header in ASCII format

    fprintf(outputFile, "BEGIN_HEADER\n");
    fprintf(outputFile, "HDR_VERSION = 1.0\n");
    fprintf(outputFile, "DATATYPE = 4D_SU3_GAUGE_3x3\n");
    fprintf(outputFile, "STORAGE_FORMAT = 1.0\n");
    fprintf(outputFile, "DIMENSION_1 = %d\n", CommonPrms::instance()->Lx());
    fprintf(outputFile, "DIMENSION_2 = %d\n", CommonPrms::instance()->Ly());
    fprintf(outputFile, "DIMENSION_3 = %d\n", CommonPrms::instance()->Lz());
    fprintf(outputFile, "DIMENSION_4 = %d\n", CommonPrms::instance()->Lt());
    fprintf(outputFile, "LINK_TRACE = %1.15g\n",link_trace);
    fprintf(outputFile, "PLAQUETTE = %1.15g\n",plaquette );
    fprintf(outputFile, "CHECKSUM = %x\n", chksum );
    fprintf(outputFile, "ENSEMBLE_ID = %s\n", "JLQCD_null");
    fprintf(outputFile, "SEQUENCE_NUMBER = %d\n", 0 );
    fprintf(outputFile, "FLOATING_POINT = IEEE64BIG\n");
    fprintf(outputFile, "END_HEADER\n");


    
  }


  int NERSCWriter::write (double *buffer, unsigned int size){
    //write data as it is (now, no compression)
    int status = EXIT_SUCCESS;
    
#ifndef BIG_ENDIAN_TYPE
    byte_swap(buffer, size);
#endif
    
    size_t res = fwrite(buffer, sizeof(double) ,size , outputFile);
    return res;
    

  }

  NERSCWriter::NERSCWriter(int tot_vol,int tot_in,int tot_ex, const Field& U):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex),
    writer_first_call(1),
    num_blocks(tot_vol/CommonPrms::instance()->Nx()){
    Staples stpl;
    plaquette = stpl.plaquette(GaugeField(U));
    link_trace = stpl.link_trace(GaugeField(U));
    chksum = checksum(GaugeField(U), 18);
  };


}
