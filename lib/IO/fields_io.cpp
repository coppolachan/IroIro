/*! 
 * @file fields_io.cpp 
 *
 * @brief Declarations of MPI safe read/write routines for fields
 *
 * Time-stamp: <2013-09-12 16:39:53 cossu>
 */

#include "fields_io.hpp"

namespace CCIO {
  QCDheader* NOheader(FILE *inputFile, fpos_t *pos){};

  GeneralReader* getReader(const std::string readerID,
			   int total_vol,
			   int Nin,
			   int Nex){
    if (readerID.compare("Binary")==0)
      return new BinaryReader(total_vol, Nin, Nex);
    if (readerID.compare("NERSC")==0)
      return new NERSCReader(total_vol, Nin, Nex);
    if (readerID.compare("JLQCDLegacy")==0)
      return new JLQCDLegacyReader(total_vol, Nin, Nex);
    if (readerID.compare("ILDG")==0)
      return new ILDGReader(total_vol, Nin, Nex);
  }

  GeneralWriter* getWriter(const std::string writerID,
			   int total_vol,
			   int Nin,
			   int Nex){
    if (writerID.compare("Binary")==0)
      return new BinaryWriter(total_vol, Nin, Nex); 
    if (writerID.compare("ILDG")==0)
      return new ILDGWriter(total_vol, Nin, Nex); 
  }

}
