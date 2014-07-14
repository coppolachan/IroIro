/*! 
 * @file NERSC_writer.cpp 
 *
 * @brief Definition of the NERSCWriter class methods
 *
 * Time-stamp: <2014-07-14 17:21:24 neo>
 */

#include "NERSC_writer.hpp"
#include "include/messages_macros.hpp"
#include "Tools/byteswap.hpp"
#include "include/errors.hpp"
#include "Measurements/GaugeM/staples.hpp"
#include <string.h>

namespace CCIO {
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
    fprintf(outputFile, "DATATYPE = 4D_SU3_GAUGE\n");
    fprintf(outputFile, "STORAGE_FORMAT = 1.0\n");
    fprintf(outputFile, "DIMENSION_1 = %d\n", CommonPrms::instance()->Lx());
    fprintf(outputFile, "DIMENSION_2 = %d\n", CommonPrms::instance()->Ly());
    fprintf(outputFile, "DIMENSION_3 = %d\n", CommonPrms::instance()->Lz());
    fprintf(outputFile, "DIMENSION_4 = %d\n", CommonPrms::instance()->Lt());
    fprintf(outputFile, "LINK_TRACE = %f\n",0.0 );
    fprintf(outputFile, "PLAQUETTE = %f\n",plaquette );
    fprintf(outputFile, "CHECKSUM = %x\n", 0 );
    fprintf(outputFile, "ENSEMBLE_ID = %s\n", "JLQCD_null");
    fprintf(outputFile, "SEQUENCE_NUMBER = %d\n", 0 );

    fprintf(outputFile, "END_HEADER");


    
  }


  int NERSCWriter::write (double *buffer, unsigned int size){
    //write data as it is (now, no compression)
    int status = EXIT_SUCCESS;
    if (writer_first_call){
      writer_first_call = 0;
 
      /* Write record header */
      header();
    }

    //write data on file
    size_t res = fwrite(buffer, sizeof(double) ,size , outputFile);
    return res;


  }

  NERSCWriter::NERSCWriter(int tot_vol,int tot_in,int tot_ex, const Field& U):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex),
    writer_first_call(1),
    num_blocks(tot_vol/CommonPrms::instance()->Nx()){
    Staples stpl;
    double plaquette = stpl.plaquette(GaugeField(U));
  };


}
