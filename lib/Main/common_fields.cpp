/*!
 * @file common_fields.cpp
 * @brief function definitions in FieldUtils
 */
#include "include/common_fields.hpp"
#include "Tools/sunMatUtils.hpp"

namespace FieldUtils{
  const Field TracelessAntihermite(const GaugeField& fce){

    using namespace SUNmat_utils;
    
    int Ndim = CommonPrms::instance()->Ndim();
    int Nvol = CommonPrms::instance()->Nvol();
    
    Field fp(fce.Format.size());
    for(int mu=0; mu< Ndim; ++mu){
      for(int site=0; site<Nvol; ++site){
	SUNmat a_h = u(fce,site,mu);
	fp.set(fce.Format.cslice(0,site,mu),anti_hermite(a_h));
      }
    }
    return fp;
  }

  GaugeField1DType DirSlice(const GaugeFieldType& F, int dir){
    return GaugeField1DType(Field(F.data[F.format.dir_slice(dir)]));
  }


}
