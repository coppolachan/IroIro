/*! 
 * @file binary_reader.cpp 
 *
 * @brief Definition of the BinaryReader class methods
 *
 * Time-stamp: <2013-09-12 12:49:47 cossu>
 */


#include "binary_reader.hpp"
#include "Tools/sunMat.hpp"
#include "include/errors.hpp"
#include "include/commonPrms.hpp"

namespace CCIO {
  void BinaryReader::set_sources(FILE *src){
    inputFile = src;
  }
  int BinaryReader::format(int gsite, int in, int ex) const{
      return in +total_internal*(ex+ total_external*gsite);
  }

  int BinaryReader::read(double *buffer, unsigned int size){
    //Read just as it is
    //order is (in, ex, sites)  
    size_t res = fread(buffer,sizeof(double),size,inputFile);
    return res;
  }

  int BinaryReader::header(){};

  int BinaryReader::check(Field& F){
    int local_volume = CommonPrms::Nvol();
    // should add a unitarity check if gauge field
    for (int mu = 0; mu < total_external; mu++){
      for (int i = 0; i < local_volume; i++){
	int local_idx = total_internal*(i +local_volume*mu);// format_G
	std::slice mat_slice = std::slice(local_idx,total_internal,1);
	SUNmat m(F[mat_slice]);

	if (!m.is_unitary()) {
	  std::ostringstream msg;
	  msg << "Failed at site = "<< i << " mu = "<< mu;
	  Errors::BaseWarning ("BinaryReader::check unitarity fail",  msg);
	  return CHECK_ERROR;
	}

      }
    }

    return CHECK_PASS;

  }; 

  BinaryReader::BinaryReader(int tot_vol,int tot_in,int tot_ex):
    total_volume(tot_vol),total_internal(tot_in),
    total_external(tot_ex){};
}


